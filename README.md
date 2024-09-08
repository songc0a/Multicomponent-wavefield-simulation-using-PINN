# Multicomponent-wavefield-simulation-using-PINN
This repository gives the codes for multicomponent wavefield simulation using Physics-informed Neural Networks.
**This repository reproduces the results of the paper "[Simulating Multicomponent Elastic Seismic Wavefield Using Deep Learning.](https://ieeexplore.ieee.org/abstract/document/10054624)"  IEEE Geoscience and Remote Sensing Letters 20, 3001105.

# Overview

We propose to use physics Informed Neural Network (PINN) to solve for multicomponent scattered wavefields for the Helmholtz equation. 
![du](https://github.com/songc0a/Multicomponent-wavefield-simulation-using-PINN/assets/31889731/7167b40c-1206-48cc-9139-3ed843b77d20)
The real part of scattered wavefields of the horizontal displacement obtained from (a) the finite-difference method, (b) from PINN. (c) Scattered wavefield difference between figures a and b.
![dv](https://github.com/songc0a/Multicomponent-wavefield-simulation-using-PINN/assets/31889731/58e15262-d3ea-48d2-b4e3-e12da122c040)
The real part of scattered wavefields of the vertical displacement obtained from (a) the finite-difference method, (b) from PINN. (c) Scattered wavefield difference between figures a and b.

# Installation of Tensorflow1

CPU usage: pip install --pre "tensorflow==1.15.*"

GPU usage: pip install --pre "tensorflow-gpu==1.15.*"

# Code explanation

helm_pinn_elastic_sx.py: Tensorflow code for solving the multicomponent scattered wavefields using PINN  
helm_pinn_elastic_sx_inversion.ipynb: Solving the multicomponent scattered wavefields and predicting corresponding P- and S-wave velocities
Elastic_sigsbee_trainingdata_generation.m: Matlab code for generating training and test data  

# Citation information

If you find our codes and publications helpful, please kindly cite the following publications.

@article{song2023simulating,
  title={Simulating Multicomponent Elastic Seismic Wavefield Using Deep Learning},
  
  author={Song, Chao and Liu, Yang and Zhao, Pengfei and Zhao, Tianshuo and Zou, Jingbo and Liu, Cai},
  
  journal={IEEE Geoscience and Remote Sensing Letters},
  
  volume={20},
  
  pages={1--5},
  
  year={2023},
  
  publisher={IEEE}
}

# contact information
If there are any problems, please contact me through my emails: chao.song@kaust.edu.sa;chaosong@jlu.edu.cn
