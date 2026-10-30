"""Audited RivGraph/CartoDEM/terrain directions for shared river profiles.

Unknown is a first-class result. Shortest-path guesses alone never establish a
direction. Topological continuity can propagate accepted evidence; cycles and
unexplained interior sinks/sources invalidate implicated links.
"""
import numpy as np
import networkx as nx
from scipy.stats import theilslopes

from .vendor.rivgraph import dir_shortest_paths_links, guess_synthetic_slope


def reach_height(reference,reach,distance,min_drop=3.):
    r=np.asarray(reach); x=np.asarray(distance);good=np.isfinite(reference[r])
    y=reference[r[good]];x=x[good]
    if len(y)<3 or x[-1]-x[0]<90:
        return 0,dict(reason='Too little longitudinal bank evidence')
    # Bound quadratic Theil-Sen work on very long channels, deterministically.
    take=np.unique(np.linspace(0,len(y)-1,min(len(y),128)).round().astype(int))
    slope,_,lower,upper=theilslopes(y[take],x[take],alpha=.9)
    drop=float(-slope*(x[-1]-x[0]))
    report=dict(drop_m=drop,slope_interval=[float(lower),float(upper)],stations=len(y))
    if abs(drop)<min_drop or lower<=0<=upper:
        return 0,report
    # Reuse the fork's attribute-slope guess explicitly with measured robust
    # CartoDEM levels, rather than its default synthetic outlet distance.
    links=dict(id=[0],conn=[[0,1]],len=[float(x[-1]-x[0])],sinuosity=[1.],guess=[[]],guess_alg=[[]])
    nodes=dict(id=[0,1],bank_level=[drop,0.])
    guess_synthetic_slope(links,nodes,node_column='bank_level',alg=100,min_slope=0,max_slope=np.inf)
    return (1 if links['guess'][0][0]==0 else -1),report


def reach_routing(network,routing,width):
    """Project normalized pixel-flow vectors onto the local centreline tangent."""
    n=len(network['station_pixels'])
    summed=np.zeros((n,2));weight=np.zeros(n);count=np.zeros(n,int)
    owner_station=np.where(network['fraction']<=.5,network['left'],network['right'])
    good=owner_station>=0
    for axis in range(2):
        np.add.at(summed[:,axis],owner_station[good],routing['vector'][good,axis])
    np.add.at(weight,owner_station[good],routing['weight'][good])
    np.add.at(count,owner_station[good],routing['count'][good])
    results=[]
    for reach in network['reaches']:
        r=np.asarray(reach);pixels=network['station_pixels'][r]
        xy=np.column_stack((pixels//width,pixels%width)).astype(float)
        tangent=np.gradient(xy,axis=0)
        norm=np.linalg.norm(tangent,axis=1)
        dot=np.sum(summed[r]*tangent,axis=1)/np.maximum(norm*weight[r],1e-12)
        usable=(count[r]>=3)&(norm>0)&(np.abs(dot)>=.15)
        positive=int(np.sum(usable&(dot>0)));negative=int(np.sum(usable&(dot<0)))
        agreement=max(positive,negative)/max(1,positive+negative)
        direction=0
        if positive+negative>=3 and agreement>=.75:
            direction=1 if positive>negative else -1
        results.append((direction,dict(forward_stations=positive,reverse_stations=negative,agreement=agreement)))
    return results


def infer_directions(network,reference,coast_station,*,routing=None,width=None,spacing=180.):
    reaches=network['reaches'];distances=network['distances'];nr=len(reaches)
    graph=nx.MultiGraph()
    for i,(r,x) in enumerate(zip(reaches,distances)):
        graph.add_edge(int(r[0]),int(r[-1]),key=i,weight=float(x[-1]))
    direction=np.zeros(nr,dtype=np.int8)
    diagnostics=[dict(reach=i,endpoints=[int(r[0]),int(r[-1])],status='unresolved',evidence=[]) for i,r in enumerate(reaches)]
    height=[]
    flow=reach_routing(network,routing,width) if routing is not None else [(0,{})]*nr
    for i,(r,x) in enumerate(zip(reaches,distances)):
        guess,detail=reach_height(reference,r,x)
        height.append(guess);diagnostics[i]['bank']=detail;diagnostics[i]['terrain_routing']=flow[i][1]
    # Ocean contacts may be interior stations: use metric distance on the full
    # station graph to identify terminal mouths without crossing land.
    sg=nx.Graph()
    for r,x in zip(reaches,distances):
        sg.add_weighted_edges_from((int(u),int(v),float(d)) for u,v,d in zip(r[:-1],r[1:],np.diff(x)))
    coast_distance=nx.multi_source_dijkstra_path_length(sg,[int(v) for v in coast_station if v in sg]) if any(v in sg for v in coast_station) else {}
    outlets=set();shortest={}
    for members in nx.connected_components(graph):
        terminals=sorted(v for v in members if graph.degree(v)==1)
        sea=[v for v in terminals if coast_distance.get(v,np.inf)<=2*spacing]
        low=[]
        if not sea:
            supported=sorted((float(reference[v]),v) for v in terminals if np.isfinite(reference[v]))
            if len(supported)>=2 and supported[1][0]-supported[0][0]>=3:
                low=[supported[0][1]]
        local_outlets=sea or low
        outlets.update(sea)
        if not local_outlets:
            continue
        ids=sorted(i for u,v,i in graph.edges(members,keys=True))
        node_ids=sorted(members)
        links=dict(id=ids,conn=[[int(reaches[i][0]),int(reaches[i][-1])] for i in ids],
                   len=[float(distances[i][-1]) for i in ids],guess=[[] for _ in ids],guess_alg=[[] for _ in ids])
        # Build incidence once; repeated list.index scans become cubic on a
        # large connected river, despite each junction having only a few links.
        incident={v:[] for v in node_ids}
        for i,(a,b) in zip(ids,links['conn']):
            incident[a].append(i);incident[b].append(i)
        nodes=dict(id=node_ids,conn=[incident[v] for v in node_ids],
                   outlets=local_outlets,inlets=[],idx=node_ids)
        dir_shortest_paths_links(links,nodes,difthresh=1.)
        for j,i in enumerate(ids):
            if links['guess'][j]:
                shortest[i]=1 if links['guess'][j][0]==reaches[i][0] else -1
    for i,r in enumerate(reaches):
        h=height[i];f=flow[i][0];s=shortest.get(i,0)
        diagnostics[i].update(bank_direction=int(h),terrain_direction=int(f),outlet_path_direction=int(s))
        mouths=[v for v in (r[0],r[-1]) if v in outlets]
        if len(mouths)==1:
            direction[i]=1 if mouths[0]==r[-1] else -1
            diagnostics[i]['evidence']=['GLO-30 ocean mouth']
            if h and h!=direction[i]:
                diagnostics[i]['bank_conflict_at_ocean']=True
        elif len(mouths)>1:
            diagnostics[i]['reason']='Both ends meet the ocean; tidal direction unresolved'
        elif h and f and h!=f:
            diagnostics[i]['reason']='CartoDEM bank trend conflicts with Terrain routing'
        elif h:
            direction[i]=h;diagnostics[i]['evidence']=['CartoDEM robust bank slope via RivGraph attribute guess']
            if f==h: diagnostics[i]['evidence'].append('Terrain routing agrees')
        elif f and f==s:
            direction[i]=f;diagnostics[i]['evidence']=['Terrain routing and RivGraph outlet-path agreement']
    # RivGraph continuity rule, implemented iteratively with evidence conflict
    # checks (its recursive setter also propagates width/parallel assumptions).
    for _ in range(nr+1):
        changed=False
        for node in sorted(graph.nodes):
            incident=sorted(i for _,_,i in graph.edges(node,keys=True))
            unknown=[i for i in incident if direction[i]==0]
            known=[i for i in incident if direction[i]!=0]
            if len(unknown)!=1 or not known: continue
            incoming=[(reaches[i][-1] if direction[i]>0 else reaches[i][0])==node for i in known]
            if len(set(incoming))!=1: continue
            i=unknown[0]
            if diagnostics[i].get('reason'): continue
            # If all known links enter, the remaining link must leave, and vice versa.
            desired=(1 if reaches[i][0]==node else -1)*(1 if incoming[0] else -1)
            if any(g and g!=desired for g in (height[i],flow[i][0])): continue
            direction[i]=desired;diagnostics[i]['evidence']=['Network continuity from accepted directions'];changed=True
        if not changed: break
    # A directed cycle has no physically downhill profile. Do not randomly flip
    # a link to force a DAG; invalidate every link internal to that SCC.
    directed=nx.DiGraph()
    for i,r in enumerate(reaches):
        if direction[i]: directed.add_edge(int(r[0] if direction[i]>0 else r[-1]),int(r[-1] if direction[i]>0 else r[0]))
    cyclic=[c for c in nx.strongly_connected_components(directed) if len(c)>1]
    bad=set()
    for c in cyclic:
        bad.update(i for i,r in enumerate(reaches) if r[0] in c and r[-1] in c)
    for i in bad:
        direction[i]=0;diagnostics[i]['reason']='Directed cycle; requires review'
    for node in sorted(graph.nodes):
        incident=sorted(i for _,_,i in graph.edges(node,keys=True))
        if len(incident)<2 or any(direction[i]==0 for i in incident): continue
        incoming=[(reaches[i][-1] if direction[i]>0 else reaches[i][0])==node for i in incident]
        if len(set(incoming))==1:
            for i in incident:
                direction[i]=0;diagnostics[i]['reason']='Unexplained interior source or sink'
    for i in range(nr):
        diagnostics[i]['direction']=int(direction[i])
        diagnostics[i]['status']='accepted' if direction[i] else 'unresolved'
    return direction,diagnostics
