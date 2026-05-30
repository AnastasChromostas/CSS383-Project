#PubChem data processing
import pubchempy as pcp
import pandas as pd

#Smiles String processing
drugs = ["Fluoxetine", "Escitalopram", "Sertraline", "Desvenlafaxine", "Duloxetine", 
         "Venlafaxine", "Amoxapine", "Amitriptyline", "Desipramine", "Trazodone", 
         "Nefazodone", "Mirtazapine", "Bupropion", "Selegiline", "Isocarboxazid", 
         "Phenelzine", "Esketamine", "Brexanolone", "Gepirone", "Zuranolone"]

#List for storing the final data
data = []

for drug in drugs:
    #Gets the compound data using the drug name 
    results = pcp.get_compounds(drug, 'name')
    
    if results:
        smiles_string = results[0].connectivity_smiles
        #Saves data as dictionaries for column separation in the chart  
        data.append({"Drug:": drug, "SMILES": smiles_string})
        #print("SMILES determined for " + drug)
    else:
        print("Data not found for " + drug) 

#Converts data into a table 
df = pd.DataFrame(data)
#Exports the table to a CSV file 
df.to_csv("antidepressant_smiles.csv", index = False)

print("SMILES strings Complete")


#ClinPGx data processing
df_anns = pd.read_csv("summary_annotations.tsv", sep="\t")
df_alleles = pd.read_csv("summary_ann_alleles.tsv", sep ="\t")

target_evidence = ["1A", "1B", "2", "3"]
#'|' = or os it finds any of the target drugs 
drug_pattern = '|'.join(drugs)
target_phenotypes = ["Toxicity", "Metabolism/PK"]

#Filters the drugs based on name, phenotype toxicty, and level of evidence 
df_anns_filtered = df_anns[
    (df_anns["Drug(s)"].str.contains(drug_pattern, case=False, na=False)) & 
    (df_anns["Phenotype Category"].isin(target_phenotypes)) & 
    (df_anns["Level of Evidence"].isin(target_evidence))].copy()

#Scoring the risk for genotypes based on annotation text
def score_genotype(row):
    text = str(row["Annotation Text"]).lower()
    #High risk
    if any(word in text for word in ["increased risk", "poor metabolizer", "higher risk", "increased severity", "increased chance"]):
        return 2
    elif "intermediate" in text or "reduced" in text:
    #Moderate risk
        return 1
    #Low/no risk
    else:
        return 0

#Uses the function to calculate the risk score for all the data in alleles 
df_alleles["Risk_Score"] = df_alleles.apply(score_genotype, axis = 1)
#Merges all the relevant data into one data frame 
df_anns_final = pd.merge(
    df_anns_filtered, 
    df_alleles[["Summary Annotation ID", "Genotype/Allele", "Annotation Text", "Risk_Score"]],
    on = "Summary Annotation ID", 
    how="inner")

final_columns = ["Drug(s)", "Gene", 
                 "Variant/Haplotypes", "Level of Evidence", 
                 "Genotype/Allele", "Risk_Score", "Annotation Text"]
df_export = df_anns_final[final_columns]

df_export.to_csv("annotations_data.csv", index=False)
print(f"Processed {len(df_export)} annotations")
#print(df_export.head())

#openFDA data processing 
import requests
import time

reaction_data = []
url = "https://api.fda.gov/drug/event.json"
skip = ["COMPLETED SUICIDE", "DRUG INEFFECTIVE", "OFF LABEL USE", 
        "TOXICITY TO VARIOUS AGENTS", "PAIN", "INTENTIONAL OVERDOSE", 
        "SUICIDE ATTEMPT", "SUICIDE IDEATION", "DRUG ABUSE", "OVERDOSE",
       "FEELING ABNORMAL", "INTENTIONAL SELF-INJURY", "PRODUCT USE IN UNAPPROVED INDICATION", 
        "DEPRESSION", "DEPRESSED MOOD", "DRUG EXPOSURE DURING PREGNANCY", "INJURY",
        "MATERNAL EXPOSURE DURING PREGNANCY", "CONDITION AGGRAVATED", "DRUG INTERACTION", "SUICIDAL IDEATION", "ANXIETY"]
for drug in drugs:
    #&count= groups by specific reaction names
    query = f'search=patient.drug.medicinalproduct.exact:"{drug.upper()}"+AND+serious:1&count=patient.reaction.reactionmeddrapt.exact&limit=100'
    
    try:
        response = requests.get(f"{url}?{query}")

        if response.status_code == 200:
            data = response.json()

            #Counts the number of distinct reactions for one drug 
            count = 0
            if "results" in data:
                for i in data["results"]:
                    reaction = i["term"].upper()

                    #We want to gather 50 (clean) ones
                    #clean = not in the skip list
                    if reaction not in skip and count <50:
                        reaction_data.append({
                            "Drug: ": drug,
                            "Reaction: ": reaction,
                            "Count": i["count"]
                        })
                        count+=1
                #print(f"Collected reactions for {drug}")
            else:
                print(f"No reaction data for {drug}")
        elif response.status_code == 404:
            print(f"No results for {drug}")
            
    except Exception as e:
        print(f"Error for {drug}: {e}")
        
    time.sleep(0.25)

df_fda = pd.DataFrame(reaction_data)
df_fda.to_csv("openfda_data.csv", index=False)

print("Complete")
print(df_fda.head(5))
