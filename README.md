# Multicomponent-wavefield-simulation-using-PINN
This repository gives the codes for multicomponent wavefield simulation using Physics-informed Neural Networks.
**This repository reproduces the results of the paper "[Simulating Multicomponent Elastic Seismic Wavefield Using Deep Learning.](https://ieeexplore.ieee.org/abstract/document/10054624)"  IEEE Geoscience and Remote Sensing Letters 20, 3001105.

# Overview

We propose to use the Fourier feature physics Informed Neural Network (PINN) to solve for multifrequency multisource scattered wavefields for the Helmholtz equation. The proposed method breaks the limitation of the numerical solver in single-frequency wavefield simulation. The structure of the Fourier feature PINN is shown below.
![FFPINN-en](https://github.com/songc0a/Fourier-feature-PINN-based-multifrequency-multisource-Helmholtz-solver/assets/31889731/35539da5-41c8-4fa5-bfb0-23fd066a3cfc)

For the velocity in (a), the scattered wavefields (5-10 Hz) from the finite difference are shown in (b), and the scattered wavefields (5-10 Hz) from the  Fourier feature PINN are shown in (c).

# Installation of Tensorflow1

CPU usage: pip install --pre "tensorflow==1.15.*"

GPU usage: pip install --pre "tensorflow-gpu==1.15.*"

# Code explanation

helm_pinn_elastic_sx.py: Tensorflow code for solving the multicomponent scattered wavefields using PINN  
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
