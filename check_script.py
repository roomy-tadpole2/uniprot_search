def check(search_key_word: str, organism: str, taxonomy_id: int):
    print(f"checking {search_key_word}, taxonomy_id={taxonomy_id}, organism={organism}...")
    result = uniprot_search(search_key_word, taxonomy_id=taxonomy_id)

    data = []
    for entry in result:
        protein_id = entry.get("primaryAccession", "N/A")
        annotation = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A")
        sequence = entry.get("sequence", "").get("value", "")
        organism = entry.get("organism", "").get("scientificName", "")
        data.append({
                "ID": protein_id,
                "Annotation": annotation,
                "Organism": organism,
                "Sequence": sequence
                })
    data = pd.DataFrame(data)
    display(data)