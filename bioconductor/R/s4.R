
# Talk about S4 classes

# S4 classes represetns complicated data structures
library(ALL)
library(GenomicRanges)

df <- data.frame(y = rnorm(10), x = rnorm(10))
line <- lm(y ~ x, data = df)

names(line) #
class(line)

# Turn R objcets to class S4
xx <- list(a = c(1,2,3), b = rnorm(3))
class(xx) <- 'lm'

isS4(xx)
isS4(ALL)

# help
class?ExpressionSet

# OLD way of create new classes
new("ExpressionSet")

# Get definition of class
getClass("ExpressionSet")

# Slots
# Where the data is. Access it with @ or slot()
# You should use the access function though =>
annotation(ALL)
# Should be documented in the help page for the class

# Name

# Extends
# OOP things


validObject(ALL)
# Check if an object is proper

##################################################
# S4 methods
##################################################
library(GenomicRanges)
# Allows you to different thing depending on the input type (shadowing)

# The method is marked as standardGeneric

showMethods('as.data.frame')
# Get all the classes for which the function is implemented
getMethod("as.data.frame", signature(x= 'GenomicRanges'))

# Each method can have its own help page ...

method?"as.data.frame,DataFrame"


showMethods("findOverlaps")
# Dispatching on both input arguments
getMethod("findOverlaps", signature(query = 'Ranges', subject = "Ranges"))

?"findOverlaps,Ranges,Ranges-method"













