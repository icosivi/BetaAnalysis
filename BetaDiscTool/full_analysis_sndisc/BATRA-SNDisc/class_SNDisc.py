import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

class SignalProbabilityModel(nn.Module):
  def __init__(self, ansatz_pmax_value):
    super().__init__()
    self.ansatz_pmax_value = ansatz_pmax_value
    self.amp_mu = nn.Parameter(torch.tensor(2.0*ansatz_pmax_value))   # Langaus centre
    self.amp_sigma = nn.Parameter(torch.tensor(3.0)) # Langaus width
    self.pow_mu = nn.Parameter(torch.tensor(2.0*ansatz_pmax_value))  # Langaus centre
    self.pow_sigma = nn.Parameter(torch.tensor(3.0)) # Langaus width
    self.invjitter_mu = nn.Parameter(torch.tensor(4.0*ansatz_pmax_value))   # Langaus centre
    self.invjitter_sigma = nn.Parameter(torch.tensor(3.0)) # Langaus width

    self.area_mu = nn.Parameter(torch.tensor(2.0*ansatz_pmax_value))   # Langaus centre
    self.area_sigma = nn.Parameter(torch.tensor(1.0)) # Langaus width

    self.noise_amp_mu = nn.Parameter(torch.tensor(0.5 * ansatz_pmax_value))
    self.noise_amp_sigma = nn.Parameter(torch.tensor(1.0))
    self.noise_pow_mu = nn.Parameter(torch.tensor(0.5 * ansatz_pmax_value))
    self.noise_pow_sigma = nn.Parameter(torch.tensor(1.0))
    self.noise_invjitter_mu = nn.Parameter(torch.tensor(0.7 * ansatz_pmax_value))
    self.noise_invjitter_sigma = nn.Parameter(torch.tensor(1.0))
    self.noise_area_mu = nn.Parameter(torch.tensor(0.5 * ansatz_pmax_value))
    self.noise_area_sigma = nn.Parameter(torch.tensor(1.0))

  def forward(self, amplitudes, pow, invjitter, area):
    amp_prob = torch.sigmoid((amplitudes - self.amp_mu) / self.amp_sigma)
    pow_prob = torch.sigmoid((pow - self.pow_mu) / self.pow_sigma)
    invjitter_prob = torch.sigmoid((invjitter - self.invjitter_mu) / self.invjitter_sigma)
    area_prob = torch.sigmoid((area - self.area_mu) / self.area_sigma)
    signal_score = amp_prob * pow_prob #* area_prob

    noise_amp_pdf = torch.exp(-0.5 * ((amplitudes - self.noise_amp_mu) / self.noise_amp_sigma) ** 2)
    noise_amp_pdf = noise_amp_pdf / (self.noise_amp_sigma * np.sqrt(2 * np.pi))
    noise_pow_pdf = torch.exp(-0.5 * ((pow - self.noise_pow_mu) / self.noise_pow_sigma) ** 2)
    noise_pow_pdf = noise_pow_pdf / (self.noise_pow_sigma * np.sqrt(2 * np.pi))
    noise_invjitter_pdf = torch.exp(-0.5 * ((invjitter - self.noise_invjitter_mu) / self.noise_invjitter_sigma) ** 2)
    noise_invjitter_pdf = noise_invjitter_pdf / (self.noise_invjitter_sigma * np.sqrt(2 * np.pi))
    noise_area_pdf = torch.exp(-0.5 * ((area - self.noise_area_mu) / self.noise_area_sigma) ** 2)
    noise_area_pdf = noise_area_pdf / (self.noise_area_sigma * np.sqrt(2 * np.pi))

    noise_score = noise_amp_pdf * noise_pow_pdf #* noise_area_pdf
    total_score = signal_score * (1 - noise_score)
    return torch.clamp(total_score, 1e-18, 1.0)

def differential_programming_SigProbModel(num_epochs, pmax, width, area, risetime, rms, ansatz_pmax_value, learning_rate):
  model = SignalProbabilityModel(ansatz_pmax_value)
  optimizer = optim.Adam(model.parameters(), lr=learning_rate)

  num_events = len(pmax)
  sig_events = np.count_nonzero(pmax >= ansatz_pmax_value)
  labels = np.concatenate([np.ones(sig_events), np.zeros(num_events - sig_events)])  # 1 = signal, 0 = noise

  max_amplitudes = torch.tensor(pmax, dtype=torch.float32)
  max_pow = torch.tensor(pmax / width, dtype=torch.float32)

  max_invjitter = torch.tensor(pmax * rms / risetime, dtype=torch.float32)

  safe_risetime = torch.tensor(risetime, dtype=torch.float32) + 1e-15
  max_area = torch.tensor(pmax * rms / safe_risetime, dtype=torch.float32) # test var, not actual area
  labels = torch.tensor(labels, dtype=torch.float32)

  for epoch in range(num_epochs):
    optimizer.zero_grad()
    scores = model(max_amplitudes, max_pow, max_invjitter, max_area)
    loss = -torch.mean(labels * torch.log(scores + 1e-12) + (1 - labels) * torch.log(1 - scores + 1e-12))  
    loss.backward()
    optimizer.step()
    if epoch % 2000 == 0:
      print(f"Epoch {epoch}, Loss: {loss.item()}")
      #for name, param in model.named_parameters():
      #  print(name, param.grad)

  return model(max_amplitudes, max_pow, max_invjitter, max_area).detach().numpy(), max_amplitudes
