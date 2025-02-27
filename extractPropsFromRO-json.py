import json
import requests
import csv

# URL of the JSON file
url = "https://raw.githubusercontent.com/oborel/obo-relations/master/ro.json"

# Fetch the JSON data
response = requests.get(url)
data = response.json()

# Extract inverse relationships
inverse_relations = {}

# Navigate to edges in the JSON structure
edges = data.get('graphs', [])[0].get('edges', [])

for edge in edges:
    subject = edge.get('sub', '')  # forward
    predicate = edge.get('pred', '')  # forward
    inverse_of = edge.get('obj', '') # inverse
    if predicate == "inverseOf":
        if inverse_of:
            inverse_relations[subject] = inverse_of

# Print results
for relation, inverse in inverse_relations.items():
    print(f"{relation}, {inverse}")

outputFile = "relations.txt"
with open(outputFile, mode='w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    for relation, inverse in inverse_relations.items():
        writer.writerow([relation, inverse])

