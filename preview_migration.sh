#!/bin/bash
# Dry-run preview: Shows what will be changed without modifying anything
# Run this from the genophi/ repository root directory

echo "=== GenoPHI Migration Dry-Run Preview ==="
echo "This will show what would be changed WITHOUT making any modifications"
echo ""

# Check if phage_modeling directory exists
if [ ! -d "phage_modeling" ]; then
    echo "ERROR: phage_modeling/ directory not found!"
    exit 1
fi

echo "Files that would be deleted:"
echo "  - phage_modeling.egg-info/"
echo "  - phage_modeling/__pycache__/"
echo "  - phage_modeling/workflows/__pycache__/"
echo ""

echo "Python files that would be updated:"
find phage_modeling/ -name "*.py" -type f | while read file; do
    count=$(grep -c "phage_modeling" "$file" 2>/dev/null || echo "0")
    if [ "$count" -gt 0 ]; then
        echo "  - $file ($count references)"
    fi
done
echo ""

echo "setup.py changes:"
grep -n "phage_modeling" setup.py 2>/dev/null || echo "  (none found)"
echo ""

echo "COPYRIGHT changes:"
grep -n "phage_modeling" COPYRIGHT 2>/dev/null || echo "  (none found)"
echo ""

echo "Directory rename:"
echo "  phage_modeling/ → genophi/"
echo ""

echo "Total Python files to update:"
find phage_modeling/ -name "*.py" -type f | wc -l
echo ""

echo "Estimated 'phage_modeling' occurrences in source code:"
grep -r "phage_modeling" phage_modeling/ --include="*.py" 2>/dev/null | wc -l
echo ""

echo "=== End of Dry-Run ==="
echo ""
echo "If this looks correct, run: ./migrate_to_genophi.sh"
