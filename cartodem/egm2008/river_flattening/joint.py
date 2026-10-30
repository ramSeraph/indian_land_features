"""Joint river profiles: one fit for all connected, supported sections.

Shared junction elevations are unknowns, not hard targets copied from raw water
pixels. Established directions constrain the first fit; missing directions are
proposed from the joint profile and checked in a second monotone fit. Those
model-derived directions are explicitly distinguished from independent evidence.
"""
import numpy as np
import networkx as nx
from scipy import sparse
from scipy.sparse.csgraph import dijkstra, connected_components

from .directions import infer_directions
from .profiles import fit_network
from .boundaries import native_values, read_contacts


def assign_stations(network, owners):
    """Map native skeleton owners to the nearest longitudinal station."""
    ids = network['ids']
    indices = np.searchsorted(ids, owners)
    safe = np.minimum(indices,len(ids)-1)
    good = (indices<len(ids)) & (ids[safe]==owners)
    stations = np.where(network['fraction'][safe]<=.5,network['left'][safe],network['right'][safe])
    good &= stations>=0
    return stations,good


def collapse_levels(n,edges,flat_reaches):
    """Equality joins share one unknown, including exactly level directed loops."""
    equality=[(r[0],v) for r in flat_reaches for v in r[1:]]
    eq=np.asarray(equality,dtype=int).reshape(-1,2)
    matrix=sparse.csr_matrix((np.ones(len(eq)),(eq[:,0],eq[:,1])),shape=(n,n))
    _,mapping=connected_components(matrix,directed=False)
    edge=np.asarray(edges,dtype=int).reshape(-1,2)
    contracted=mapping[edge]
    count=int(mapping.max())+1
    matrix=sparse.csr_matrix((np.ones(len(edge)),(contracted[:,0],contracted[:,1])),shape=(count,count))
    _,cycle=connected_components(matrix,directed=True,connection='strong')
    mapping=cycle[mapping]
    contracted=mapping[edge]
    contracted=np.unique(contracted[contracted[:,0]!=contracted[:,1]],axis=0)
    return mapping,contracted


def complete_components(network,tables):
    """Never publish a fitted fragment beside an unresolved part of its river."""
    rejected=set();iterations=[]
    for _ in range(len(network['fitted'])+1):
        values,groups=native_values(network);bad=set();pairs=0
        for table in tables:
            owners,_=read_contacts(table)
            if not len(owners):continue
            index=np.searchsorted(network['ids'],owners)
            if np.any(index>=len(values)) or np.any(network['ids'][index]!=owners):
                raise ValueError('Contact owner is absent from joint network')
            active=np.isfinite(values[index]);g=groups[index]
            unresolved=(active[:,0]!=active[:,1])|(active.all(axis=1)&(g[:,0]!=g[:,1]))
            pairs+=int(unresolved.sum())
            bad.update(int(v) for v in g[unresolved][active[unresolved]] if v>=0)
        if not bad:
            return dict(method='Shared junction variables; complete connected treatment',
                        raw_height_pins=False,boundary_pairs=0,rejected_profile_groups=sorted(rejected),
                        iterations=iterations)
        iterations.append(dict(unresolved_contacts=pairs,rejected_groups=sorted(bad)))
        rejected.update(bad)
        drop=np.isin(network['component'],list(bad));network['fitted'][drop]=np.nan
        a,b=network['left'],network['right'];good=(a>=0)&(b>=0)
        network['native_supported'][good]&=~(drop[a[good]]|drop[b[good]])
        edges=network['supported_edges'];network['supported_edges']=edges[~(drop[edges[:,0]]|drop[edges[:,1]])]
    raise ValueError('Connected-treatment validation did not converge')


def fit_joint(native,network,*,routing=None,width=None,spacing=180.,quantile=.2,
              stiffness=2.,max_support_gap=1200.,contacts=()):
    n=len(network['station_pixels']);os,good=assign_stations(network,native['bank_owners'])
    os=os[good];oh=native['bank_heights'][good]
    cs,good=assign_stations(network,native['ocean_owners']);cs=np.unique(cs[good])
    count=np.bincount(os,minlength=n);order=np.argsort(os,kind='stable');starts=np.r_[0,np.cumsum(count)]
    reference=np.full(n,np.nan)
    for i in np.flatnonzero(count>=4):reference[i]=np.quantile(oh[order[starts[i]:starts[i+1]]],quantile)
    seeds,diagnostics=infer_directions(network,reference,cs,routing=routing,width=width,spacing=spacing)
    graph=nx.Graph()
    for r,x in zip(network['reaches'],network['distances']):
        graph.add_weighted_edges_from((int(a),int(b),float(d)) for a,b,d in zip(r[:-1],r[1:],np.diff(x)))
    fitted=np.full(n,np.nan,np.float32);labels=np.full(n,-1,np.int32)
    reports=[];final_edges=[];applied=set()
    for group,members_set in enumerate(sorted(nx.connected_components(graph),key=min)):
        members=np.array(sorted(members_set));labels[members]=group;observed=np.isfinite(reference[members])
        ri=[i for i,r in enumerate(network['reaches']) if r[0] in members_set]
        report=dict(component=group,stations=len(members),reaches=len(ri),supported_stations=int(observed.sum()),status='unchanged')
        reports.append(report)
        if len(members)<3 or observed.sum()<3 or observed.mean()<.6:
            report['reason']='Insufficient bank support for connected river';continue
        adjacency=nx.to_scipy_sparse_array(graph.subgraph(members_set),nodelist=members.tolist(),weight='weight',format='csr')
        adjacency=sparse.csr_matrix(adjacency);adjacency.indices=adjacency.indices.astype('int32');adjacency.indptr=adjacency.indptr.astype('int32')
        gap=dijkstra(adjacency,directed=False,indices=np.flatnonzero(observed),min_only=True)
        report['maximum_support_gap_m']=float(gap.max())
        if gap.max()>max_support_gap:
            report['reason']='Connected river has an unsupported section';continue
        reaches=[np.searchsorted(members,network['reaches'][i]).tolist() for i in ri]
        distances=[network['distances'][i] for i in ri]
        seed_edges=[]
        for i,r in zip(ri,reaches):
            if seeds[i]:
                rr=r if seeds[i]>0 else list(reversed(r));seed_edges.extend(zip(rr[:-1],rr[1:]))
        sample=np.isin(os,members);localobs=np.searchsorted(members,os[sample]);z=oh[sample]
        caps=np.where(observed,reference[members],np.inf)
        weights=1/np.bincount(localobs,minlength=len(members))[localobs]
        try:
            initial,_=fit_network(len(members),seed_edges,localobs,z,upper_bounds=caps,
                                  quantile=quantile,stiffness=stiffness,weights=weights,
                                  unoriented_reaches=reaches,unoriented_distances=distances)
            edges=[];flat=[];directions={}
            for i,r in zip(ri,reaches):
                sign=int(seeds[i])
                if not sign:
                    difference=float(initial[r[0]])-float(initial[r[-1]])
                    if abs(difference)<=.001:
                        flat.append(r);directions[i]=0;continue
                    sign=1 if difference>0 else -1
                directions[i]=sign;rr=r if sign>0 else list(reversed(r));edges.extend(zip(rr[:-1],rr[1:]))
            mapping,contracted=collapse_levels(len(members),edges,flat)
            size=int(mapping.max())+1;upper=np.full(size,np.inf)
            np.minimum.at(upper,mapping,caps)
            final,detail=fit_network(size,contracted,mapping[localobs],z,upper_bounds=upper,
                                     quantile=quantile,stiffness=stiffness,weights=weights,
                                     unoriented_reaches=[mapping[r].tolist() for r in reaches],
                                     unoriented_distances=distances)
            values=final[mapping]
            change=np.abs(values[observed]-reference[members][observed])
            report.update(p95_station_adjustment_m=float(np.percentile(change,95)),maximum_station_adjustment_m=float(change.max()))
            if np.percentile(change,95)>15 or change.max()>40:
                raise ValueError('Joint profile requires excessive adjustment of supported bank levels')
            # Check the full profile of every reach after shared-level expansion,
            # rather than relying on endpoint sign alone.
            for i,r in zip(ri,reaches):
                delta=np.diff(values[r]);direction=directions[i]
                if (direction and np.any(delta*direction>0)) or (not direction and np.any(delta!=0)):
                    raise ValueError('Expanded connected profile violates its downstream/level constraints')
            fitted[members]=values;applied.update(ri)
            # Use original station indices for audit, including flat links.
            for i,r in zip(ri,reaches):
                rr=r if directions[i]>=0 else list(reversed(r))
                final_edges.extend((int(members[a]),int(members[b])) for a,b in zip(rr[:-1],rr[1:]))
                diagnostics[i]['joint_direction']=directions[i]
                diagnostics[i]['direction_source']='independent evidence' if seeds[i] else 'joint profile hypothesis'
            report.update(status='fitted',fit=detail,joint_direction_reaches=sum(not seeds[i] for i in ri),
                          level_reaches=len(flat),shared_level_unknowns=size)
        except ValueError as exc:report['reason']=str(exc)
    left,right=network['left'],network['right'];supported=(left>=0)&(right>=0)
    supported[supported]&=np.isfinite(fitted[left[supported]])&np.isfinite(fitted[right[supported]])
    network.update(fitted=fitted,component=labels,native_supported=supported,
                   supported_edges=np.asarray(final_edges,dtype=int).reshape(-1,2),bank_reference=reference,bank_count=count,
                   direction_diagnostics=diagnostics)
    boundary=complete_components(network,list(contacts));network['boundary_acceptance']=boundary
    rejected=set(boundary['rejected_profile_groups'])
    for report in reports:
        if report['component'] in rejected:
            report.update(status='unchanged',reason='Connected river still contains unmodelled native contacts')
    for i,d in enumerate(diagnostics):
        d['profile_applied']=i in applied and labels[d['endpoints'][0]] not in rejected
    return network,reports
