import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
from sentence_transformers import SentenceTransformer, util
import random



class ContextualTaxonomyDataset(Dataset):
    def __init__(self, triplets):
        self.triplets = triplets  # (anchor, positive, negative) triplets
    def __len__(self):
        return len(self.triplets)
    def __getitem__(self, idx):
        return self.triplets[idx]


class ContextualTMN(nn.Module):
    def __init__(self, embedding_model='paraphrase-MiniLM-L6-v2'):
        super(ContextualTMN, self).__init__()
        self.encoder = SentenceTransformer(embedding_model)  # pre-trained model
        self.loss_fn = nn.TripletMarginLoss(margin=1.0)  # triplet Loss
    def forward(self, anchor, positive, negative):
        # Convert to tensor with requires_grad=True
        anchor_emb = self.encoder.encode(anchor, convert_to_tensor=True).requires_grad_()
        positive_emb = self.encoder.encode(positive, convert_to_tensor=True).requires_grad_()
        negative_emb = self.encoder.encode(negative, convert_to_tensor=True).requires_grad_()
        loss = self.loss_fn(anchor_emb, positive_emb, negative_emb)
        return loss

# sample Triplet Data (anchor, positive, negative)
triplets = [
    ("Drosophila parasitoidOf Apis mellifera", 
     "Drosophila (Fungi) parasitoidOf Apis mellifera", 
     "Drosophila melanogaster (Fruit Fly) modelOrganism"),
    ("Candida infects humans", 
     "Candida albicans (Pathogenic Yeast) infects humans", 
     "Candida tropicalis (Non-pathogenic) lives in water")
]

device = "cuda" if torch.cuda.is_available() else "cpu" #device
model = ContextualTMN().to(device)

optimizer = optim.Adam(model.parameters(), lr=0.001) #optimizer

num_epochs = 10
for epoch in range(num_epochs):                     #training loops
    total_loss = 0
    for anchor, positive, negative in triplets:
        optimizer.zero_grad()  # clear previous gradients
        loss = model(anchor, positive, negative)  # forward pass
        loss.backward()  # backpropagate
        optimizer.step()  # update weights
        total_loss += loss.item()
    print(f"Epoch {epoch+1}/{num_epochs}, Loss: {total_loss/len(triplets)}")


def find_best_match(query_name, taxonomy_list, model, device='cpu'):
    model.to(device)
    query_emb = model.encode(query_name, convert_to_tensor=True, device=device)
    taxonomy_embs = model.encode(taxonomy_list, convert_to_tensor=True, device=device)
    query_emb = query_emb.to(device)
    taxonomy_embs = taxonomy_embs.to(device)
    similarities = util.pytorch_cos_sim(query_emb, taxonomy_embs) #cosine similarity
    best_match_idx = similarities.argmax().item() #best match
    best_match = taxonomy_list[best_match_idx]
    best_score = similarities[0][best_match_idx].item()
    return best_match, best_score

#prelim db with context
taxonomy_db = [
    "Drosophila melanogaster (Fruit Fly) modelOrganism",
    "Drosophila (Fungi) parasitoidOf Apis mellifera",
    "Candida albicans (Pathogenic Yeast) infects humans",
    "Candida tropicalis (Non-pathogenic) lives in water"
]

# initialize model and move to CUDA if available
device = 'cuda' if torch.cuda.is_available() else 'cpu'
model = SentenceTransformer('paraphrase-MiniLM-L6-v2').to(device)
query_name = "Drosophila parasitoidOf Apis mellifera"
print("db:")
print(taxonomy_db)
print("query:")
print(query_name)
best_match, score = find_best_match(query_name, taxonomy_db, model, device)

print(f"best Match: {best_match}, similarity score: {score:.4f}")


