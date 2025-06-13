import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import struct
import numpy as np

batch_size = 128
filename = "training.data"
input_dim = 4  # Each element is a 4-vector
phi_dim = 128
rho_dim = 64
gamma_dim = 32
eon_size = 500


class DeepElementSelector(nn.Module):
    def __init__(self, input_dim, phi_dim, rho_dim, gamma_dim):
        super(DeepElementSelector, self).__init__()
        # Phi network: simplified with only LayerNorm
        self.phi = nn.Sequential(
            nn.Linear(input_dim, phi_dim),
            nn.ReLU(),
            nn.LayerNorm(phi_dim),
            nn.Linear(phi_dim, phi_dim),
            nn.ReLU(),
            nn.LayerNorm(phi_dim),
            nn.Linear(phi_dim, phi_dim // 2),
        )

        # Rho network: processes aggregated representation
        self.rho = nn.Sequential(
            nn.Linear(phi_dim // 2, rho_dim),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(rho_dim, rho_dim),
            nn.ReLU(),
            nn.LayerNorm(rho_dim),
            nn.Linear(rho_dim, rho_dim // 2),
        )

        # Gamma network: simplified with only LayerNorm
        self.gamma = nn.Sequential(
            nn.Linear(input_dim, gamma_dim),
            nn.ReLU(),
            nn.LayerNorm(gamma_dim),
            nn.Linear(gamma_dim, gamma_dim),
            nn.ReLU(),
            nn.LayerNorm(gamma_dim),
            nn.Linear(gamma_dim, rho_dim // 2),  # Same dim as rho output
        )

    def forward(self, x):
        batch_size, num_elements, input_dim = x.shape
        x_flat = x.view(-1, input_dim)  # [batch_size * num_elements, input_dim]
        phi_out = self.phi(x_flat).view(batch_size, num_elements, -1)
        aggregated = torch.sum(phi_out, dim=1)  # [batch_size, phi_dim // 2]
        global_context = self.rho(aggregated).unsqueeze(
            1
        )  # [batch_size, 1, rho_dim // 2]
        element_transformed = self.gamma(x_flat).view(batch_size, num_elements, -1)
        scores = torch.sum(element_transformed * global_context, dim=2)
        selection_probs = F.softmax(scores, dim=1)  # [batch_size, num_elements]
        return {
            "scores": scores,
            "selection_probs": selection_probs,
        }


def read_batch_from_binary_file(file_handle):
    # Read set size (4 bytes, int32)
    set_size = struct.unpack("<i", file_handle.read(4))[0]
    num_doubles = set_size * input_dim * batch_size
    doubles = []

    for i in range(num_doubles):
        double_bytes = file.read(8)  # 8 bytes per double
        if len(double_bytes) < 8:
            break
        double_value = float(struct.unpack("d", double_bytes)[0])
        doubles.append(double_value)

    # Convert to numpy array and reshape
    shaped_data = np.array(doubles, dtype=np.float32).reshape(
        (batch_size, set_size, input_dim)
    )
    return (set_size, torch.from_numpy(shaped_data))


# Initialize the model and optimizer (back to float32)

model = DeepElementSelector(input_dim, phi_dim, rho_dim, gamma_dim)
optimizer = torch.optim.Adam(model.parameters(), lr=0.0005)

# Training loop with binary file reading
with open(filename, "rb") as file:
    print("training network to select first compton scatter\n")
    epoch = 0
    accuracy = 0
    while True:
        # Read batch from binary file
        set_size, shaped_data = read_batch_from_binary_file(file)

        # Target is always the first element (index 0) in each set
        target_indices = torch.zeros(batch_size, dtype=torch.long)

        optimizer.zero_grad()

        # Forward pass
        outputs = model(shaped_data)
        dot_products = outputs["scores"]
        selection_probs = outputs["selection_probs"]

        # Loss function: use dot_products for cross_entropy
        loss = F.cross_entropy(dot_products, target_indices)

        # Backward pass
        loss.backward()
        optimizer.step()

        with torch.no_grad():
            predicted_indices = torch.argmax(dot_products, dim=1)
            accuracy += (predicted_indices == target_indices).float().mean()
        if epoch % eon_size == 0:
            accuracy /= eon_size
            # Calculate accuracy for monitoring

            print(f"Eon {int(epoch / eon_size)}:")
            print(f"  Accuracy: {accuracy.item():.4f}")
            print()
            accuracy = 0

        epoch += 1
