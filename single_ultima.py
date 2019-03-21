#importing custom module for analysis
import convert as cvt
import scipy.misc as sm
import numpy as np

#example to create SPEI image for one image

#create 100 x 100 matrix of ones
img_a = np.ones((100,100))
#save image
sm.imsave('image_a.png', img_a) #result should look completely black

#setup image for finding circle
img_b = sm.imread('image_a.png', flatten=True)
img_c = cvt.convert(img_b)
#save image
sm.imsave('image_c.png', img_c) #result should have a white square with large
                                #boarder

#encircling image
ultima = cvt.enc_circ(img_c)
#save final image
sm.imsave('ultima.png', ultima) #result should be same white square with smaller
                                #boarder than img_c

#



#
