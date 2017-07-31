import copy
"""
GlobalFunction is a function needed to use multiprocess with a function of an object.
It takes as input: the object, the name of the function and the input of the function.
It makes a copy of the object, take the function and finally launch the function.
"""
def GlobalFunction(inputAll):
    (objectClass,functionName, inputFunction) = inputAll
    objectCopy = copy.copy(objectClass)
    function = getattr(objectCopy,functionName)
    if type(inputFunction) is tuple: output = function(*inputFunction)
    else:   output = function(inputFunction)
    del objectCopy,inputAll,objectClass,functionName, inputFunction,function
    return output

