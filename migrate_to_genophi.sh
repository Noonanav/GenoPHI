#!/bin/bash
# Migration script: phage_modeling → genophi
# Run this from the genophi/ repository root directory

set -e  # Exit on error

echo "=== GenoPHI Package Renaming Script ==="
echo ""

# Pre-flight checks
echo "Running pre-flight checks..."

# Check if phage_modeling directory exists
if [ ! -d "phage_modeling" ]; then
    echo "ERROR: phage_modeling/ directory not found!"
    echo "Are you in the correct directory?"
    exit 1
fi

# Check if genophi directory already exists
if [ -d "genophi" ]; then
    echo "ERROR: genophi/ directory already exists!"
    echo "Please remove it first or this is already migrated."
    exit 1
fi

# Check if setup.py exists
if [ ! -f "setup.py" ]; then
    echo "ERROR: setup.py not found!"
    echo "Are you in the repository root?"
    exit 1
fi

# Warn about git status
if [ -d ".git" ]; then
    if ! git diff-index --quiet HEAD -- 2>/dev/null; then
        echo "⚠️  WARNING: You have uncommitted changes in git"
        echo "It's recommended to commit or stash changes before migration"
        read -p "Continue anyway? (y/n) " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Aborted. Please commit your changes first."
            exit 1
        fi
    fi
fi

echo "✓ Pre-flight checks passed"
echo ""

# Final confirmation
read -p "This will rename phage_modeling to genophi. Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 1
fi

echo ""
echo "Step 1: Cleaning up generated files..."
rm -rf phage_modeling.egg-info/
rm -rf phage_modeling/__pycache__/
rm -rf phage_modeling/workflows/__pycache__/
echo "✓ Removed .egg-info and __pycache__ directories"

echo ""
echo "Step 2: Updating imports in ALL Python files in package..."
# Update all .py files in phage_modeling (including workflows and core modules)
find phage_modeling/ -name "*.py" -type f -exec sed -i.bak 's/from phage_modeling\./from genophi./g' {} +
find phage_modeling/ -name "*.py" -type f -exec sed -i.bak 's/import phage_modeling\./import genophi./g' {} +
find phage_modeling/ -name "*.py" -type f -exec sed -i.bak 's/from phage_modeling import/from genophi import/g' {} +
find phage_modeling/ -name "*.py" -type f -exec sed -i.bak 's/import phage_modeling$/import genophi/g' {} +
echo "✓ Updated all Python imports in phage_modeling/"

echo ""
echo "Step 3: Updating setup.py..."
sed -i.bak "s/name='phage_modeling'/name='genophi'/g" setup.py
sed -i.bak 's/phage_modeling\.workflows/genophi.workflows/g' setup.py
echo "✓ Updated setup.py"

echo ""
echo "Step 4: Updating COPYRIGHT..."
sed -i.bak 's/phage_modeling/genophi/g' COPYRIGHT
echo "✓ Updated COPYRIGHT"

echo ""
echo "Step 5: Renaming package directory..."
mv phage_modeling genophi
echo "✓ Renamed phage_modeling/ → genophi/"

echo ""
echo "Step 6: Cleaning up backup files..."
find . -name "*.bak" -type f -delete
echo "✓ Removed .bak files"

echo ""
echo "=== Migration Complete! ==="
echo ""

# Verification
echo "Verifying migration..."
if [ -d "genophi" ] && [ ! -d "phage_modeling" ]; then
    echo "✓ Directory renamed successfully"
else
    echo "⚠️  Warning: Unexpected directory state"
fi

# Check for remaining references (excluding README, manuscript_scripts, and .git)
remaining=$(grep -r "phage_modeling" genophi/ --exclude-dir=".git" --exclude-dir="manuscript_scripts" --exclude="*.ipynb" --exclude="README.md" 2>/dev/null | grep -v "Binary file" | wc -l)
if [ "$remaining" -eq 0 ]; then
    echo "✓ No remaining 'phage_modeling' references in source code"
else
    echo "⚠️  Found $remaining remaining references to 'phage_modeling'"
    echo "    (This might be expected in comments or manuscript files)"
fi

echo ""
echo "Next steps:"
echo "1. Review changes: git diff"
echo "2. Test import: python -c 'import genophi; print(genophi.__version__)'"
echo "3. Update README.md (manual - references to phage_modeling)"
echo "4. Reinstall package: pip install -e ."
echo "5. Run a simple test workflow to verify"
echo "6. Commit changes: git add -A && git commit -m 'Rename package to genophi'"
echo "7. Rename GitHub repo: Settings → Repository name → genophi"
echo "8. Update remote: git remote set-url origin git@github.com:Noonanav/genophi.git"
echo ""

echo ""
echo "=== Migration Complete! ==="
echo ""
echo "Next steps:"
echo "1. Review changes: git diff"
echo "2. Test imports: python -c 'import genophi'"
echo "3. Update README.md (references to phage_modeling)"
echo "4. Rename GitHub repo: Settings → Repository name → genophi"
echo "5. Update remote URL: git remote set-url origin git@github.com:Noonanav/genophi.git"
echo "6. Reinstall package: pip install -e ."
echo ""
