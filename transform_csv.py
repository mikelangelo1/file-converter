import csv
import sys

def split_csv_columns(input_file, output_file, group_size=16000):
    # Read the input CSV file
    with open(input_file, 'r') as f:
        reader = csv.reader(f)
        data = [row[0] for row in reader]  # Assuming single column

    # Calculate the number of columns needed
    num_columns = (len(data) + group_size - 1) // group_size

    # Prepare the output data
    output_data = []
    for i in range(min(group_size, len(data))):
        row = []
        for j in range(num_columns):
            index = j * group_size + i
            if index < len(data):
                row.append(data[index])
            else:
                row.append('')  # Fill empty cells with blank
        output_data.append(row)

    # Write the output CSV file
    if sys.version_info[0] < 3:
        # Python 2
        with open(output_file, 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(output_data)
    else:
        # Python 3
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(output_data)


# Example usage
input_file = './mutant_variants_copy.csv'
output_file = './mutant_variants_transformed.csv'
split_csv_columns(input_file, output_file)

# print(f"CSV file has been processed. Output saved to {output_file}")