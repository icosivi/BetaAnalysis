for file in dirOfExampleCards/Run4_*.txt; do
    echo "Running analysis on $file ..."
    python3 analyseInPY3.py "$file"

    # Check exit status
    if [ $? -ne 0 ]; then
        echo "Analysis failed for $file, skipping to next."
    fi
done

echo "All analyses attempted!"
