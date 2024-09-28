import pandas as pd

def calculate_coverage_stats(file_path):
    # Read the coverage file with the correct delimiter and stripping whitespace
    coverage_df = pd.read_csv(file_path, sep='\s+', header=None, names=['Scaffold', 'Position', 'Coverage'])

    # Print the first few rows to inspect the data
    print("First few rows of the data:")
    print(coverage_df.head())

    # Check the data types
    print("\nData types:")
    print(coverage_df.dtypes)

    # Convert Coverage column to numeric, coercing errors to NaN
    coverage_df['Coverage'] = pd.to_numeric(coverage_df['Coverage'], errors='coerce')

    # Drop rows where Coverage is NaN
    coverage_df = coverage_df.dropna(subset=['Coverage'])

    # Calculate statistics
    max_coverage = coverage_df['Coverage'].max()
    min_coverage = coverage_df['Coverage'].min()
    mean_coverage = coverage_df['Coverage'].mean()
    quartiles = coverage_df['Coverage'].quantile([0.25, 0.5, 0.75, 0.8, 0.9, 0.95])

    # Print the results
    print(f"\nMax Coverage: {max_coverage}")
    print(f"Min Coverage: {min_coverage}")
    print(f"Mean Coverage: {mean_coverage}")
    print("Quartiles:")
    print(quartiles)

    return max_coverage, min_coverage, mean_coverage, quartiles

if __name__ == "__main__":
    import sys
    file_path = sys.argv[1]
    calculate_coverage_stats(file_path)
