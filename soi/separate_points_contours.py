#!/usr/bin/env python3
"""
Separate GeoJSONL file into points (circles) and contours.
Points are identified as circles with similar radius and elevation not a multiple of 10.
"""

import json
import sys
from shapely.geometry import shape, Point
import numpy as np


def is_circle(geometry, tolerance=0.03):
    """
    Check if a geometry is a circle by verifying all points are equidistant from centroid.
    
    Args:
        geometry: Shapely geometry object
        tolerance: Relative tolerance for radius variation (default 15%)
    
    Returns:
        tuple: (is_circle: bool, radius: float or None)
    """
    if geometry.geom_type != 'LineString':
        return False, None
    
    coords = list(geometry.coords)
    if len(coords) < 8:  # Too few points to be a meaningful circle
        return False, None
    
    # Calculate centroid
    centroid = geometry.centroid
    
    # Calculate distances from centroid to all points
    distances = []
    for coord in coords:
        point = Point(coord)
        distance = centroid.distance(point)
        distances.append(distance)
    
    distances = np.array(distances)
    mean_radius = np.mean(distances)
    
    if mean_radius == 0:
        return False, None
    
    # Calculate coefficient of variation (std/mean)
    std_dev = np.std(distances)
    cv = std_dev / mean_radius

    #print(cv, mean_radius, file=sys.stderr)
    
    # If coefficient of variation is low, it's likely a circle
    if cv < tolerance:
        return True, mean_radius
    
    return False, None


def separate_geojsonl(input_file, points_file, contours_file):
    """
    Separate GeoJSONL file into points and contours.
    
    Args:
        input_file: Input GeoJSONL file path
        points_file: Output file for points (circles)
        contours_file: Output file for contours
    """
    points_count = 0
    contours_count = 0
    truncation_count = 0
    
    with open(input_file, 'r') as infile, \
         open(points_file, 'w') as points_out, \
         open(contours_file, 'w') as contours_out:
        
        for line_num, line in enumerate(infile, 1):
            if line_num % 10000 == 0:
                print(f"Processed {line_num} features... (Points: {points_count}, Contours: {contours_count})", file=sys.stderr)
            
            try:
                feature = json.loads(line.strip())
                
                # Get elevation value
                elevation = float(feature['properties'].get('VALUE', 0))
                
                feature['properties']['VALUE'] = elevation
                
                # Convert to shapely geometry
                geom = shape(feature['geometry'])
                
                # Check if it's a circle
                is_circ, radius = is_circle(geom)
                
                if is_circ and radius < 0.0001:
                    # print('Found circle:', is_point, 'with radius:', radius, 'and elevation:', elevation, file=sys.stderr)
                    # Convert to point at centroid
                    centroid = geom.centroid
                    point_feature = {
                        'type': 'Feature',
                        'properties': {
                            **feature['properties'],
                            'radius': radius,
                        },
                        'geometry': {
                            'type': 'Point',
                            'coordinates': [centroid.x, centroid.y]
                        },
                        'id': feature.get('id')
                    }
                    points_out.write(json.dumps(point_feature) + '\n')
                    points_count += 1
                else:
                    # Keep as contour
                    contours_out.write(line)
                    contours_count += 1
                    
            except Exception as e:
                print(f"Error processing line {line_num}: {e}", file=sys.stderr)
                # Write to contours by default on error
                contours_out.write(line)
                contours_count += 1
    
    print(f"\nComplete!", file=sys.stderr)
    print(f"Points (circles): {points_count}", file=sys.stderr)
    print(f"Contours: {contours_count}", file=sys.stderr)
    print(f"Elevation values truncated by > 0.01: {truncation_count}", file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python separate_points_contours.py <input.geojsonl> <points.geojsonl> <contours.geojsonl>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    points_file = sys.argv[2]
    contours_file = sys.argv[3]
    
    print(f"Processing {input_file}...", file=sys.stderr)
    separate_geojsonl(input_file, points_file, contours_file)
