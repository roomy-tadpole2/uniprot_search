import requests

def uniprot_search(query, taxonomy_id, size=100):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"{query} AND taxonomy_id:{taxonomy_id}",
        "format": "json",
        "size": size  # Limit
    }

    response = requests.get(base_url, params=params)

    if response.status_code == 200:
        data = response.json()
        results = data.get("results", [])
        
        if not results:
            return {"error": "No protein data found for the given query."}
        
        return results
    else:
        return print({
            "error": f"Failed to fetch data. HTTP Status Code: {response.status_code}",
            "response_text": response.text
        })

def parse_data(results, organism_:str, keywords_list:list):
    
    banwords_list = ['transferase','demethylase']
    parsed_data = []

    for entry in results:
        
        protein_id = entry.get("primaryAccession", "N/A")
        annotation = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A")
        sequence = entry.get("sequence", "").get("value", "")
        organism = entry.get("organism", "").get("scientificName", "")

        parse_flag = True
        if (organism_ not in organism): parse_flag=False
        for keyword in keywords_list:
            if (keyword.lower() not in annotation.lower()): parse_flag=False
        for banword in banwords_list:
            if (banword in annotation): parse_flag=False
        if (parse_flag):
            parsed_data.append({
                "ID": protein_id,
                "Annotation": annotation,
                "Organism": organism,
                "Sequence": sequence
            })

    return parsed_data