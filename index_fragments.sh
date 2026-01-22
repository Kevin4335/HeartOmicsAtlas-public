#!/bin/bash
# index_fragments.sh
# Install tabix and index ATAC fragment files

echo "=========================================="
echo "Indexing ATAC Fragment Files"
echo "=========================================="

# Check if tabix is installed
if ! command -v tabix &> /dev/null; then
    echo "Step 1: Installing tabix (htslib)..."
    sudo apt-get update
    sudo apt-get install -y tabix
    if [ $? -ne 0 ]; then
        echo "Error: Failed to install tabix. Trying alternative method..."
        # Try installing via conda or other method
        echo "Please install tabix manually: sudo apt-get install tabix"
        exit 1
    fi
    echo "✓ tabix installed successfully"
else
    echo "✓ tabix is already installed"
fi

# Check tabix version
echo ""
echo "tabix version:"
tabix 2>&1 | head -1

# Index fragment files
echo ""
echo "Step 2: Indexing fragment files..."
BASE_DIR="/home/ubuntu/HeartOmics/resources-NEW/multi_omics_data"

for dir in 3655 3675 7668 3688; do
    frag_file="$BASE_DIR/$dir/atac_fragments.tsv.gz"
    index_file="$BASE_DIR/$dir/atac_fragments.tsv.gz.tbi"
    
    if [ -f "$frag_file" ]; then
        if [ -f "$index_file" ]; then
            echo "  ✓ $dir/atac_fragments.tsv.gz.tbi already exists"
        else
            echo "  Indexing $dir/atac_fragments.tsv.gz..."
            cd "$BASE_DIR/$dir"
            tabix -p bed atac_fragments.tsv.gz
            if [ $? -eq 0 ] && [ -f "atac_fragments.tsv.gz.tbi" ]; then
                echo "    ✓ Successfully indexed"
            else
                echo "    ✗ Failed to index"
            fi
        fi
    else
        echo "  ✗ $frag_file not found"
    fi
done

echo ""
echo "=========================================="
echo "Indexing complete!"
echo "=========================================="
echo ""
echo "Verifying index files..."
for dir in 3655 3675 7668 3688; do
    index_file="$BASE_DIR/$dir/atac_fragments.tsv.gz.tbi"
    if [ -f "$index_file" ]; then
        size=$(ls -lh "$index_file" | awk '{print $5}')
        echo "  ✓ $dir/atac_fragments.tsv.gz.tbi exists ($size)"
    else
        echo "  ✗ $dir/atac_fragments.tsv.gz.tbi missing"
    fi
done
