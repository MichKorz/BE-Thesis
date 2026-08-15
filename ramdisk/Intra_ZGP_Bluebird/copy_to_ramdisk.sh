#!/bin/bash

TARGET_DIR="/mnt/ramdisk"
SCRIPT_NAME=$(basename "$0")

mkdir -p "$TARGET_DIR"

# Iteracja po wszystkich plikach/folderach w obecnym katalogu
for item in *; do
    # Kopiowanie wszystkiego z wyjątkiem samego skryptu
    if [ "$item" != "$SCRIPT_NAME" ]; then
        cp -r "$item" "$TARGET_DIR/"
    fi
done

echo "Kopiowanie zakończone."
