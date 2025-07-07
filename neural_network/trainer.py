import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import struct
import numpy as np
import sys

batch_size = 128
input_dim = 4  # Each element is a 4-vector
phi_dim = 128
rho_dim = 64
gamma_dim = 32
eon_size = 500


class DeepElementSelector(nn.Module):
    def __init__(self, input_dim, phi_dim, rho_dim, gamma_dim):
        super(DeepElementSelector, self).__init__()
        self.phi = nn.Sequential(
            nn.Linear(input_dim, phi_dim),
            nn.ReLU(),
            nn.LayerNorm(phi_dim),
            nn.Linear(phi_dim, phi_dim),
            nn.ReLU(),
            nn.LayerNorm(phi_dim),
            nn.Linear(phi_dim, phi_dim // 2),
        )

        self.rho = nn.Sequential(
            nn.Linear(phi_dim // 2, rho_dim),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(rho_dim, rho_dim),
            nn.ReLU(),
            nn.LayerNorm(rho_dim),
            nn.Linear(rho_dim, rho_dim // 2),
        )

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
        global_context = self.rho(aggregated).unsqueeze(1)
        element_transformed = self.gamma(x_flat).view(batch_size, num_elements, -1)
        scores = torch.sum(element_transformed * global_context, dim=2)
        return scores


def read_batch_from_binary_file(file_handle):
    set_size_bytes = file_handle.read(4)
    if len(set_size_bytes) != 4:
        return (None, None)
    set_size = struct.unpack("<i", set_size_bytes)[0]
    num_doubles = set_size * input_dim * batch_size
    doubles = []

    for i in range(num_doubles):
        double_bytes = file.read(8)  # 8 bytes per double
        if len(double_bytes) < 8:
            break
        double_value = float(struct.unpack("d", double_bytes)[0])
        doubles.append(double_value)

    shaped_data = np.array(doubles, dtype=np.float32).reshape(
        (batch_size, set_size, input_dim)
    )
    return (set_size, torch.from_numpy(shaped_data))


if len(sys.argv) != 3:
    print(str(len(sys.argv)) + " parameters, expected 4, usage:")
    print("python3 trainer.py [training data] [output location]")
    sys.exit()
model = DeepElementSelector(input_dim, phi_dim, rho_dim, gamma_dim)
optimizer = torch.optim.Adam(model.parameters(), lr=0.0005)

with open(sys.argv[1], "rb") as file:
    print("Training network to select first compton scatter...\n")
    epoch = 0
    accuracy = 0
    set_size, shaped_data = read_batch_from_binary_file(file)
    while set_size is not None:
        target_indices = torch.zeros(batch_size, dtype=torch.long)

        optimizer.zero_grad()

        scores = model(shaped_data)

        loss = F.cross_entropy(scores, target_indices)

        loss.backward()
        optimizer.step()

        set_size, shaped_data = read_batch_from_binary_file(file)
        with torch.no_grad():
            predicted_indices = torch.argmax(scores, dim=1)
            selection_probs = F.softmax(scores, dim=1)
            accuracy += (predicted_indices == target_indices).float().mean()
        if epoch % eon_size == 0:
            accuracy /= eon_size

            print(f"Eon {int(epoch / eon_size)} | Accuracy: {accuracy.item():.4f}")
            accuracy = 0

        epoch += 1
print("All data used, saving traced neural network...")
model.eval()

# Create example input with appropriate dimensions
# Using a reasonable set_size for tracing (e.g., 50 elements)
example_set_size = 50
example_input = torch.randn(1, example_set_size, input_dim)

# Trace the model
traced_model = torch.jit.optimize_for_inference(torch.jit.trace(model, example_input))

# Save the traced model
traced_model.save(sys.argv[2])

print("Model saved as: 'chooser.pt'")
print(f"Model expects input shape: [batch_size, num_elements, {input_dim}]")
