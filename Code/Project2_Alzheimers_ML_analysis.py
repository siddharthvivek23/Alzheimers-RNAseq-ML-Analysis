import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, accuracy_score


# -------------------------------
# Load RNA-seq data
# -------------------------------

df = pd.read_csv(
    "Data/GSE53697_RNAseq_AD.txt",
    sep="\t"
)


# -------------------------------
# Prepare gene expression matrix
# -------------------------------

# Extract RPKM expression values
expression = df.iloc[:, 19:]

# Add gene names
expression.index = df["GeneSymbol"]

# Transpose:
# rows = samples
# columns = genes
X = expression.T


# -------------------------------
# Create labels
# -------------------------------

# Control = 0 (8 samples)
# Alzheimer's = 1 (9 samples)

y = [0]*8 + [1]*9


# -------------------------------
# Feature selection
# -------------------------------

# Find genes with highest variation

gene_variance = X.var(axis=0)

top_variable_genes = gene_variance.sort_values(
    ascending=False
).head(1000).index


# Keep top 1000 genes

X_selected = X[top_variable_genes]


# -------------------------------
# Train/Test split
# -------------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X_selected,
    y,
    test_size=0.2,
    random_state=42,
    stratify=y
)


# -------------------------------
# Logistic Regression Model
# -------------------------------

model = LogisticRegression(
    max_iter=1000
)

model.fit(
    X_train,
    y_train
)


# -------------------------------
# Model Evaluation
# -------------------------------

y_pred = model.predict(X_test)

cm = confusion_matrix(
    y_test,
    y_pred
)

accuracy = accuracy_score(
    y_test,
    y_pred
)

print("Confusion Matrix:")
print(cm)

print("Accuracy:")
print(accuracy)


# -------------------------------
# Extract AI Biomarker Weights
# -------------------------------

coefficients = model.coef_[0]


gene_importance = pd.DataFrame({

    "Gene": X_selected.columns,
    "Weight": coefficients

})


# Rank genes

top_genes = gene_importance.sort_values(
    by="Weight",
    ascending=False
)


print("\nTop AI Biomarkers:")
print(top_genes.head(20))


# Save biomarker list

top_genes.to_csv(
    "Final_AI_Predicted_Alzheimers_Biomarkers.csv",
    index=False
)



# -------------------------------
# Alzheimer's Gene Validation
# -------------------------------

known_AD_genes = [

    "APOE",
    "CLU",
    "APP",
    "TREM2",
    "PSEN1",
    "PSEN2",
    "BIN1",
    "ABCA7"

]


AD_validation = top_genes[
    top_genes["Gene"].isin(known_AD_genes)
]


print("\nKnown AD Genes Found:")
print(AD_validation)


# Save validation results

AD_validation.to_csv(
    "Known_AD_Genes_Identified_by_AI.csv",
    index=False
)



# -------------------------------
# Alzheimer's Gene Validation Plot
# -------------------------------

plt.figure(figsize=(8,5))


plt.bar(
    AD_validation["Gene"],
    AD_validation["Weight"]
)


plt.xlabel(
    "Alzheimer's-associated Gene"
)

plt.ylabel(
    "Logistic Regression Weight"
)

plt.title(
    "Known Alzheimer's Genes Identified by AI Model"
)


plt.axhline(0)


plt.tight_layout()

plt.show()


# -------------------------------
# GO Enrichment Analysis
# -------------------------------

from gprofiler import GProfiler


# Create gProfiler object

gp = GProfiler(
    return_dataframe=True
)


# Use top 20 AI biomarker genes

biomarker_list = top_genes.head(20)["Gene"].tolist()


print("\nAI Biomarker Genes:")
print(biomarker_list)



# Run GO enrichment

go_results = gp.profile(
    organism="hsapiens",
    query=biomarker_list
)


# Keep only Biological Process results

go_BP = go_results[
    go_results["source"] == "GO:BP"
]


print("\nTop GO Biological Processes:")
print(
    go_BP[
        [
            "name",
            "p_value"
        ]
    ].head(10)
)



# Save GO results

go_results.to_csv(
    "AI_Biomarker_GO_Enrichment.csv",
    index=False
)


# -------------------------------
# Save ML Summary Results
# -------------------------------

with open("ML_Model_Summary.txt", "w") as f:

    f.write("Alzheimer's RNA-seq Logistic Regression Model\n")
    f.write("--------------------------------------------\n\n")

    f.write("Accuracy:\n")
    f.write(str(accuracy))
    f.write("\n\n")

    f.write("Confusion Matrix:\n")
    f.write(str(cm))
    f.write("\n\n")

    f.write("Top AI Biomarkers:\n")
    f.write(top_genes.head(20).to_string())

    f.write("\n\nKnown Alzheimer's Genes Identified:\n")
    f.write(AD_validation.to_string())


# -------------------------------
# Save Trained ML Model
# -------------------------------

import joblib

joblib.dump(
    model,
    "Alzheimers_LogisticRegression_Model.pkl"
)



# -------------------------------
# Top AI Biomarker Heatmap
# -------------------------------

import seaborn as sns

# Select top 20 AI biomarker genes
heatmap_genes = top_genes.head(20)["Gene"]

# Extract expression values
heatmap_data = X[heatmap_genes]

# Create heatmap

plt.figure(figsize=(10,6))

sns.heatmap(
    heatmap_data,
    cmap="viridis",
    center=0
)

plt.title(
    "Expression Heatmap of Top AI Alzheimer's Biomarkers"
)

plt.xlabel(
    "AI Biomarker Genes"
)

plt.ylabel(
    "Samples"
)

plt.tight_layout()

plt.savefig(
    "AI_Biomarker_Heatmap.png",
    dpi=300
)

plt.show()



