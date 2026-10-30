#!/bin/bash

# This script converts all TIFF files in a source directory to
# Cloud-Optimized GeoTIFFs (COGs) in a destination directory.

# Set the source and destination directories
SOURCE_DIR="data/raw/v3_r1"
DEST_DIR="data/cogs"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Find all .tif files in the source directory and process them
for TIF_FILE in "$SOURCE_DIR"/*.tif; do
  # Get the base name of the file
  BASENAME=$(basename "$TIF_FILE")

  # Define the output file path
  OUTPUT_FILE="$DEST_DIR/$BASENAME"

  if [[ -f "$OUTPUT_FILE" ]]; then
    echo "COG already exists for $BASENAME, skipping conversion."
    continue
  fi

  # Convert the TIFF to a Cloud-Optimized GeoTIFF
  echo "Converting $TIF_FILE to $OUTPUT_FILE..."
  gdal_translate "$TIF_FILE" "$OUTPUT_FILE" -of COG -b 1 -a_nodata -255 -a_srs EPSG:4326
done

echo "All TIFF files have been converted to COGs."
