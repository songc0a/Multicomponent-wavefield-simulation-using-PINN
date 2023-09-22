import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from scipy.interpolate import griddata
import time
import cmath

np.random.seed(1234)
tf.set_random_seed(1234)


def fwd_gradients(Y, x):
    dummy = tf.ones_like(Y)
    G = tf.gradients(Y, x, grad_ys=dummy, colocate_gradients_with_ops=True)[0]
    Y_x = tf.gradients(G, dummy, colocate_gradients_with_ops=True)[0]
    return Y_x

fre = 4.0
PI = 3.1415926
niter = 50000
misfit = []
misfit1 = []
class PhysicsInformedNN:
    # Initialize the class
    def __init__(self, x, z, sx, alpha, alpha0, beta, beta0, dxxu0, dzzu0, dxzu0, dxxv0, dzzv0, dxzv0, layers):
        
        X = np.concatenate([x, z, sx], 1)
        self.iter=0
        self.start_time=0

        self.lb = X.min(0)
        self.ub = X.max(0)
                
        self.X = X
        
        self.x = X[:,0:1]
        self.z = X[:,1:2]
        self.sx = X[:,2:3]

        self.alpha = alpha
        self.alpha0 = alpha0

        self.beta = beta
        self.beta0 = beta0

        self.dxxu0 = dxxu0
        self.dzzu0 = dzzu0
        self.dxzu0 = dxzu0

        self.dxxv0 = dxxv0
        self.dzzv0 = dzzv0
        self.dxzv0 = dxzv0        

        self.layers = layers
        
        # Initialize NN
        self.weights, self.biases = self.initialize_NN(layers)  
        self.saver = tf.train.Saver() ##saver initialization

        # tf placeholders 
        self.sess = tf.Session(config=tf.ConfigProto(allow_soft_placement=True,
                                                     log_device_placement=True))
        
        self.x_tf = tf.placeholder(tf.float32, shape=[None, self.x.shape[1]])
        self.z_tf = tf.placeholder(tf.float32, shape=[None, self.z.shape[1]])
        self.sx_tf = tf.placeholder(tf.float32, shape=[None, self.sx.shape[1]])

        self.du_pred, self.dv_pred, self.fu_loss, self.fv_loss = self.net_NS(self.x_tf, self.z_tf, self.sx_tf)

        # loss function we define
        self.loss = tf.reduce_sum(tf.square(tf.abs(self.fu_loss))) + tf.reduce_sum(tf.square(tf.abs(self.fv_loss)))
        
        # optimizer used by default (in original paper)
        self.optimizer = tf.contrib.opt.ScipyOptimizerInterface(self.loss, 
                                                                method = 'L-BFGS-B', 
                                                                options = {'maxiter': 50000,
                                                                           'maxfun': 50000,
                                                                           'maxcor': 50,
                                                                           'maxls': 50,
                                                                           'ftol' : 1.0 * np.finfo(float).eps})

        self.global_step = tf.Variable(0, trainable=False)
        starter_learning_rate = 0.001
        self.learning_rate = tf.train.exponential_decay(starter_learning_rate, self.global_step, 5000, 0.9, staircase=False)

        self.optimizer_Adam = tf.train.AdamOptimizer(learning_rate=self.learning_rate)
        self.train_op_Adam = self.optimizer_Adam.minimize(self.loss)

        init = tf.global_variables_initializer()
        self.sess.run(init)
#        self.saver.restore(self.sess,'./checkpoint_dir/Mymodel_ff06_elastic_sx') # model save

    def initialize_NN(self, layers):        
        weights = []
        biases = []
        num_layers = len(layers) 
        for l in range(0,num_layers-1):
            W = self.xavier_init(size=[layers[l], layers[l+1]])
            b = tf.Variable(tf.zeros([1,layers[l+1]], dtype=tf.float32)+0.0, dtype=tf.float32)
            weights.append(W)
            biases.append(b)        
        return weights, biases
        
    def xavier_init(self, size):
        in_dim = size[0]
        out_dim = size[1]        
        xavier_stddev = np.sqrt(1./in_dim)
        return tf.Variable(tf.truncated_normal([in_dim, out_dim], stddev=xavier_stddev), dtype=tf.float32)

        # Evaluates the forward pass
    def forward_pass(self, H):
        num_layers = len(self.layers)

        X = H
        H = 2.0 * (X - self.lb) / (self.ub - self.lb) - 1.0

        for l in range(0, num_layers - 2):
            W = self.weights[l]
            b = self.biases[l]
            H = tf.atan(tf.add(tf.matmul(H, W), b))

        # Merge the outputs by concatenation

        W = self.weights[-1]
        b = self.biases[-1]
        H = tf.add(tf.matmul(H, W), b)

        return H

    # Forward pass for u
    def net_u(self, x, z, sx):
        u = self.forward_pass(tf.concat([x, z, sx], 1))
        return u

    def net_NS(self, x, z, sx):

        omega = tf.complex(2.0*PI*fre, 0.0)
        alpha = self.alpha
        alpha0 = self.alpha0
        d_alpha = alpha - alpha0
        
        beta = self.beta
        beta0 = self.beta0
        d_beta = beta - beta0

        dxxu0 = self.dxxu0
        dzzu0 = self.dzzu0
        dxzu0 = self.dxzu0

        dxxv0 = self.dxxv0
        dzzv0 = self.dzzv0
        dxzv0 = self.dxzv0

        u_and_v = self.net_u(x,z,sx)

        du_real = u_and_v[:,0:1]
        du_imag = u_and_v[:,1:2]
        dv_real = u_and_v[:,2:3]
        dv_imag = u_and_v[:,3:4]

        du = tf.complex(du_real, du_imag)

        dudx = fwd_gradients(du, x)
        dudz = fwd_gradients(du, z)

        du_xx = fwd_gradients(dudx, x)
        du_xz = fwd_gradients(dudx, z)
        du_zz = fwd_gradients(dudz, z)

        dv = tf.complex(dv_real, dv_imag)

        dvdx = fwd_gradients(dv, x)
        dvdz = fwd_gradients(dv, z)

        dv_xx = fwd_gradients(dvdx, x)
        dv_xz = fwd_gradients(dvdx, z)
        dv_zz = fwd_gradients(dvdz, z)
        
        fu_loss = omega*omega*du + alpha*du_xx + beta*du_zz + (alpha-beta)*dv_xz + d_alpha*dxxu0 + d_beta*dzzu0 + (d_alpha-d_beta)*dxzv0
        fv_loss = omega*omega*dv + alpha*dv_zz + beta*dv_xx + (alpha-beta)*du_xz + d_alpha*dzzv0 + d_beta*dxxv0 + (d_alpha-d_beta)*dxzu0

        return du, dv, fu_loss, fv_loss
    
    def callback(self, loss):
        #print('Loss: %.3e' % (loss))
        misfit1.append(loss)
        self.iter=self.iter+1
        if self.iter % 10 == 0:
                elapsed = time.time() - self.start_time
                print('It: %d, LBFGS Loss: %.3e,Time: %.2f' %
                      (self.iter, loss, elapsed))
                self.start_time = time.time()
      
    def train(self, nIter): 

        tf_dict = {self.x_tf: self.x, self.z_tf: self.z, self.sx_tf: self.sx}
 
        start_time = time.time()
        for it in range(nIter):
            self.sess.run(self.train_op_Adam, tf_dict)
            loss_value = self.sess.run(self.loss, tf_dict)
            misfit.append(loss_value)         
            # Print
            if it % 10 == 0:
                elapsed = time.time() - start_time
                loss_value = self.sess.run(self.loss, tf_dict)
                #misfit.append(loss_value)
                print('It: %d, Loss: %.3e,Time: %.2f' % 
                      (it, loss_value, elapsed))
                start_time = time.time()

            
        self.optimizer.minimize(self.sess,
                                feed_dict = tf_dict,
                                fetches = [self.loss],
                                loss_callback = self.callback)

        self.saver.save(self.sess,'./checkpoint_dir/Mymodel_ff_elastic_sx_5hz') # model save            
    
    def predict(self, x_star, z_star, sx_star):
        
        tf_dict = {self.x_tf: x_star, self.z_tf: z_star, self.sx_tf: sx_star}       
 
        du = self.sess.run(self.du_pred, tf_dict)
        dv = self.sess.run(self.dv_pred, tf_dict)

        return du, dv
        
        
if __name__ == "__main__": 

    layers = [3, 128, 128, 64, 64, 32, 32, 32, 32, 4]
    #layers = [2, 20, 20, 20, 20, 20, 20, 20, 20, 2] # neurons in each layer (used in original paper)
    
    # Load Test Data
    data = scipy.io.loadmat('Elastic_sigsbee_4Hz_sx_F10.mat')

    x_star = data['x_star'] 
    z_star = data['z_star'] 
    sx_star = data['sx_star']

    x_train = data['x_train']
    z_train = data['z_train']
    sx_train = data['sx_train']

    alpha_train = data['alpha_train']
    alpha0_train = data['alpha0_train']

    beta_train = data['beta_train']
    beta0_train = data['beta0_train']

    dxxu0_train = data['dxxu0_train']
    dzzu0_train = data['dzzu0_train']
    dxzu0_train = data['dxzu0_train']

    dxxv0_train = data['dxxv0_train']
    dzzv0_train = data['dzzv0_train']
    dxzv0_train = data['dxzv0_train']

    # Training
    model = PhysicsInformedNN(x_train, z_train, sx_train, alpha_train, alpha0_train, beta_train, beta0_train, dxxu0_train, dzzu0_train, dxzu0_train, dxxv0_train, dzzv0_train, dxzv0_train, layers)
    model.train(niter)

    # Prediction
    scipy.io.savemat('misfit_fre_elastic_%dhz_sx.mat' % (fre), {'misfit': misfit})
    scipy.io.savemat('misfit1_fre_elastic_%dhz_sx.mat' % (fre), {'misfit1': misfit1})

    du_pred, dv_pred = model.predict(x_star, z_star, sx_star)
    scipy.io.savemat('du_pred_atan_%dhz_fre_elastic_sx.mat' % (fre), {'du_pred': du_pred})
    scipy.io.savemat('dv_pred_atan_%dhz_fre_elastic_sx.mat' % (fre), {'dv_pred': dv_pred})


  
