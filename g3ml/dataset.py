from torch.utils.data import Dataset
from torchvision import transforms


# Custom dataset class
class ScDataset(Dataset):
    def __init__(self, data, transform=None):
        self.data = data
        self.transform = transform

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        return self.data[index]
