import csv

input_file = "sigma_gammap_jpsi.agr"
output_file = "sigma_gammap_jpsi.csv"

datasets = []
current_data = []
with open(input_file, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("@"):
            if current_data:
                datasets.append(current_data)
                current_data = []
            continue
        try:
            parts = line.split()
            nums = [float(p) for p in parts]
            if len(nums) == 3:
                current_data.append(nums)
        except ValueError:
            continue

if current_data:
    datasets.append(current_data)

with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["dataset", "W", "Sigma", "Error"])
    for i, data in enumerate(datasets):
        for row in data:
            writer.writerow([i] + row)

print(f"Salvo: {output_file}")
# Note: The output CSV will have columns: dataset, W, Sigma, Error