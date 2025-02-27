import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from sentence_transformers import SentenceTransformer, util

#load Taxonomy Metadata (Lineage & Interactions)
class TaxonomyMetadata:
    def __init__(self, csv_file):
        self.metadata = pd.read_csv(csv_file)  # load lineage & interactions
        self.metadata.fillna("", inplace=True)  # fill missing values with empty string

    def get_context(self, organism):
        """Retrieve lineage & interactions for a given organism."""
        match = self.metadata[self.metadata["organism"].str.lower() == organism.lower()]
        if match.empty:
            return ""
        lineage = match.iloc[0]["lineage"]
        interactions = match.iloc[0]["possible_interactions"]
        return f"Lineage: {lineage}, Possible Interactions: {interactions}"

#define Dataset Class for Context-Aware Triplets
class TaxonomyDataset(Dataset):
    def __init__(self, triplets_csv, lineage_csv):
        self.data = pd.read_csv(triplets_csv)  # load triplets
        self.metadata = TaxonomyMetadata(lineage_csv)  # load taxonomy metadata

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        anchor, relation, target = self.data.iloc[idx, 0], self.data.iloc[idx, 1], self.data.iloc[idx, 2]

        # retrieve lineage-based context
        anchor_context = self.metadata.get_context(anchor)
        target_context = self.metadata.get_context(target)

        # construct natural language descriptions
        anchor_desc = f"{anchor} ({anchor_context})"
        positive_desc = f"{anchor} {relation} {target} ({target_context})"
        negative_desc = f"{anchor} unrelated to {target} (random example)"

        return anchor_desc, positive_desc, negative_desc

#define Context-Aware TMN Model
class ContextualTMN(nn.Module):
    def __init__(self, embedding_model='paraphrase-MiniLM-L6-v2'):
        super(ContextualTMN, self).__init__()
        self.encoder = SentenceTransformer(embedding_model)
        self.loss_fn = nn.TripletMarginLoss(margin=1.0)  

    def forward(self, anchor, positive, negative):
        embeddings = self.encoder.encode(anchor + positive + negative, convert_to_tensor=True)
        anchor_emb, positive_emb, negative_emb = embeddings[:len(anchor)], embeddings[len(anchor):2*len(anchor)], embeddings[2*len(anchor):]
        return self.loss_fn(anchor_emb, positive_emb, negative_emb)

#training Loop
def train_model(triplet_csv, lineage_csv, num_epochs=10, batch_size=16, lr=0.001, device="cuda" if torch.cuda.is_available() else "cpu"):
    dataset = TaxonomyDataset(triplet_csv, lineage_csv)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    model = ContextualTMN().to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    for epoch in range(num_epochs):
        total_loss = 0
        for anchor, positive, negative in dataloader:
            optimizer.zero_grad()
            loss = model(anchor, positive, negative)  # forward
            loss.backward()  # backpropagation
            optimizer.step()  # update weights
            total_loss += loss.item()
        print(f"Epoch {epoch+1}/{num_epochs}, Loss: {total_loss/len(dataloader):.4f}")

#run Training
train_model("triplets.csv", "lineage.csv", num_epochs=10, batch_size=16)

