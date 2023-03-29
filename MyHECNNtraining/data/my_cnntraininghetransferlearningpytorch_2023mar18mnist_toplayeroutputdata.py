# -*- coding: utf-8 -*-
"""My_CNNtrainingHEtransferlearningPyTorch_2023Mar18MNIST_TopLayerOutputData.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1eHfLWYtrAmh-0aFvDE41KZOCTU4GFFwY
"""

import csv
import torch
from torch import nn
from torch.utils.data import DataLoader
from torchvision import models, datasets
import torchvision.transforms as trnsfrms
from torchvision.transforms import ToTensor, Resize, Lambda

trnsfrms = trnsfrms.Compose([Resize(224), ToTensor(),  Lambda(lambda x: x.repeat(3, 1, 1) ) ])  # Grayscale Images like MNIST and USPS
#trnsfrms = trnsfrms.Compose([Resize(224), ToTensor(), ])                                       # Color Images like CIFAR10
#trnsfrms = trnsfrms.Compose([ ToTensor(), ]) 

# Download training data from open datasets.FashionMNIST.MNIST.USPS  / CIFAR10
training_data = datasets.MNIST(
    root="data",
    train=True,
    download=True,
    transform= trnsfrms
)

# Download test data from open datasets.FashionMNIST.MNIST
testing_data = datasets.MNIST(
    root="data",
    train=False,
    download=True,
    transform= trnsfrms
)

batch_size = 512

# Create data loaders.
train_dataloader = DataLoader(training_data, batch_size=batch_size)
test_dataloader = DataLoader(testing_data, batch_size=batch_size)

# Get cpu or gpu device for training.
device = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
print(f"Using {device} device")


#### ConvNet as fixed feature extractor ####
model = models.regnet_x_400mf(pretrained=True)  #regnet_y_400mf(pretrained=True) 
for param in model.parameters():
    param.requires_grad = False
# Parameters of newly constructed modules have requires_grad=True by default
model.fc = nn.Sequential(*list(model.fc.children())[:-1])
model = model.to(device)
print(model)
print(model.fc)
print(type(model.fc))


# Save the raw dataset: USPS MNIST CIFAR10
train_dataset = []
size = len(train_dataloader.dataset)
num_batches = len(train_dataloader)
model.eval()
with torch.no_grad():
    for X, y in train_dataloader:
        X, y = X.to(device), y.to(device)
        pred = model(X)
        print(f"Shape of model(X) [N, C, H, W]: {pred.shape}")
        print(f"Shape of y: {y.shape} {y.dtype}")
        pred = torch.reshape(pred, (pred.shape[0], -1) )
        y = torch.reshape(y, (y.shape[0], -1) )
        print(f"Shape of model(X): {pred.shape} {pred.dtype}")
        print(f"Shape of y: {y.shape} {y.dtype}")

        train_dataset += torch.cat( (y, pred ), 1)
        print(f"Shape of train_dataset: {len(train_dataset)}, {len(train_dataset[0])}")

print("train_dataset :" + str(len(train_dataset)) + ",\t" + str(len(train_dataset[0])) )
print(type(train_dataset))
with open('REGNET_X_400MF_MNIST_TRAININGfirst8192.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    for i in range(8192):
      writer.writerow(train_dataset[i].detach().cpu().numpy())
csvfile.close()



# Save the raw dataset: USPS MNIST CIFAR10
test_dataset = []
size = len(test_dataloader.dataset)
num_batches = len(test_dataloader)
model.eval()
with torch.no_grad():
    for X, y in test_dataloader:
        X, y = X.to(device), y.to(device)
        pred = model(X)
        print(f"Shape of model(X) [N, C, H, W]: {pred.shape}")
        print(f"Shape of y: {y.shape} {y.dtype}")
        pred = torch.reshape(pred, (pred.shape[0], -1) )
        y = torch.reshape(y, (y.shape[0], -1) )
        print(f"Shape of model(X): {pred.shape} {pred.dtype}")
        print(f"Shape of y: {y.shape} {y.dtype}")

        test_dataset += torch.cat( (y, pred ), 1)
        print(f"Shape of test_dataset: {len(test_dataset)}, {len(test_dataset[0])}")

print("test_dataset :" + str(len(test_dataset)) + ",\t" + str(len(test_dataset[0])) )
print(type(test_dataset))
with open('REGNET_X_400MF_MNIST_TESTING.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    for i in range(len(test_dataset)):
      writer.writerow(test_dataset[i].detach().cpu().numpy())
csvfile.close()








'''
model = models.regnet_x_400mf(pretrained=True)
#print(model)
print('++++++++regnet_x_400mf++++++++++++++++++++++++++++++++++')
print(model.fc) 
#Linear(in_features=400,

model = models.regnet_y_400mf(pretrained=True)
#print(model)
print('+++++++++regnet_y_400mf+++++++++++++++++++++++++++++++++')
print(model.fc) 
#(fc): Linear(in_features=440,

model = models.regnet_x_800mf(pretrained=True)
#print(model)
print('++++++++regnet_x_800mf++++++++++++++++++++++++++++++++++')
print(model.fc) 
#Linear(in_features=672,

model = models.regnet_y_800mf(pretrained=True)
#print(model)
print('+++++++++regnet_y_800mf++++++++++++++++++++++++++++++++++')
print(model.fc) 
#(fc): Linear(in_features=784,

model = models.regnet_y_1_6gf(pretrained=True)
#print(model)
print('+++++++++regnet_y_1_6gf++++++++++++++++++++++++++++++++++')
print(model.fc) 
#Linear(in_features=888,

model = models.regnet_x_1_6gf(pretrained=True)
#print(model)
print('++++++++regnet_x_1_6gf++++++++++++++++++++++++++++++++++')
print(model.fc) 
#Linear(in_features=912,

model = models.regnet_x_3_2gf(pretrained=True)
print(model)
print('+++++++regnet_x_3_2gf+++++++++++++++++++++++++++++++++++')
print(model.fc) 
#Linear(in_features=1008,
'''