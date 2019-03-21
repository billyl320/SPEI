#importing custom module for analysis
import convert as cvt

#desired directories
#note that each class should be separated into different directories.
#however, for the fucntion to work, multiple directories should be specified.
#thus, an empty folder is utilized for this task
#the empty folder is called "none"

#----------------------------------------------------------
#EXAMPLE class

Class = ["CLASS_DIR", "none"]

#name of .txt file
name = 'class.txt'

#collecting EIs
cvt.BinaryHistTXT(name, Class)



#name of .txt file for traditional image histograms
name2 = 'class_old.txt'

#getting traditional image histogram 
cvt.BinaryHistTXT_old(name2, Class)



#
