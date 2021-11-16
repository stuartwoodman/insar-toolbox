#!/usr/bin/env bash
set -euo pipefail

INSAR_DB_NAME=${INSAR_DB_NAME:-grace}
INSAR_DB_HOST=${INSAR_DB_HOST:-localhost}
INSAR_DB_PORT=${INSAR_DB_PORT:-5432}
INSAR_DB_USER=${INSAR_DB_USER:-vgl}
# INSAR_DB_PASS must be passed in

INSAR_DB_TABLE=${INSAR_DB_TABLE:-insar_tiles}

MERGED="/tmp/merged_insar_tiles.shp"

# Merge all shapefiles into one, keeping the original layer name as the
# "tile_name" field.
ogrmerge.py -o "$MERGED" -overwrite_ds -single -src_layer_field_name tile_name -nln "$INSAR_DB_TABLE" "$@"

# Load the merged tiles into the database.
# NB this will fail if the table already exists.
ogr2ogr -f PostgreSQL PG:"dbname='$INSAR_DB_NAME' host='$INSAR_DB_HOST' port='$INSAR_DB_PORT' user='$INSAR_DB_USER' password='$INSAR_DB_PASS'" -nln "$INSAR_DB_TABLE" "$MERGED"

