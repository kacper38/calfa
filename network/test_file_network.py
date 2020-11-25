import numpy as np
import tensorflow as tf
import keras
import os
from sklearn.model_selection import train_test_split
#import sklearn

data_dir   = '/home/kacper/PycharmProjects/calfa/datasets/'
train_dir  = '/home/kacper/PycharmProjects/calfa/training/'
#%%

dataq3r =np.load(data_dir+'q3_pcacc_regular_ss.npy', allow_pickle=True)
dataq3o =np.load(data_dir+'q3_pcacc_onehot_ss.npy', allow_pickle=True)
dataq8r =np.load(data_dir+'q8_pcacc_regular_ss.npy', allow_pickle=True)
dataq8o =np.load(data_dir+'q8_pcacc_onehot_ss.npy', allow_pickle=True)
#%%
def add_padding(tab_of_vectors, final_len):
    output = []
    max_len = final_len#len(max(tab_of_vectors, key=len))
    for vector in tab_of_vectors:
        dim_diff = max_len - vector.shape[0]
        output.append(np.pad(vector, (0, dim_diff) , 'constant', constant_values=(0.0,0.0)))
    return  np.array(output)

#%%
#format is data[0] trainable,  data[1] is labels

tar = tf.keras.preprocessing.sequence.pad_sequences(dataq3r[0], padding='post', dtype='float32')
fea = tf.keras.preprocessing.sequence.pad_sequences(dataq3r[1], padding='post', dtype='int32')


X_train, X_test, y_train, y_test = train_test_split(tar, fea,
                                                    test_size=0.25, random_state=42)
batch_size=len(X_train)

#%%
max_leny = max([ len(max(x, key=len)) for x in [y_train, y_test]])
y_train = add_padding(y_train, max_leny)
y_test = add_padding(y_test, max_leny)

max_lenx = max([ len(max(x, key=len)) for x in [X_train, X_test]])
#X_train = add_padding(X_train, max_leny)
#X_test = add_padding(X_test, max_leny)
X_train_padded = tf.keras.preprocessing.sequence.pad_sequences(X_train, padding='post',dtype='float32')
X_test_padded = tf.keras.preprocessing.sequence.pad_sequences(X_test, padding='post',dtype='float32')

#%%

model = tf.keras.models.Sequential()

model.add(tf.keras.layers.Dense(3, input_dim=3852, activation=tf.nn.relu))
model.add(tf.keras.layers.Dense(1283, activation=tf.nn.sigmoid))
#%%
model.compile(loss='mean_squared_error',
              optimizer='adam',
              metrics=['binary_accuracy'])

model.summary()
#%%
history = model.fit(X_train_padded, y_train, epochs=500, verbose=2)
#%%


lc = history.history
#%%



inputs = tf.keras.Input(shape=(3852))
x = tf.keras.layers.Dense(3, activation=tf.nn.relu)(inputs)
outputs = tf.keras.layers.Dense(1283, activation=tf.nn.softmax)(x)

modelt = tf.keras.Model(inputs=inputs, outputs=outputs)
#%%
modelt.compile(optimizer='Adam', loss='mse' ,metrics=['mae','mse','accuracy'])
print(modelt.summary())
#%%
model.fit(X_train, y_train,epochs=500, verbose=2)
#%%
print(modelt.predict(X_train))
print(modelt.evaluate(X_test, y_test, verbose=2))
#%%
print(modelt.predict(X_test))
#%%
class TemporalSoftmax(keras.layers.Layer):
    def call(self, inputs, mask=None):
        broadcast_float_mask = tf.expand_dims(tf.cast(mask, "float32"), -1)
        inputs_exp = tf.exp(inputs) * broadcast_float_mask
        inputs_sum = tf.reduce_sum(inputs * broadcast_float_mask, axis=1, keepdims=True)
        return inputs_exp / inputs_sum

inputs1 = keras.Input(shape=(3852))
x1 = keras.layers.Embedding(input_dim=3852, output_dim=1283, mask_zero=True)(inputs1)
x1 = keras.layers.Dense(1)(x1)
#lstm = tf.keras.layers.LSTM(3)
outs = TemporalSoftmax()(x1)


#outputs1 = keras.layers.LSTM(1283)(x1)

modelj = keras.Model(inputs1, outputs1)
#modelj.compile(optimizer='Adam', loss='mse' ,metrics=['mae','mse','accuracy'])
#%%
modelj.fit(X_train, y_train)
