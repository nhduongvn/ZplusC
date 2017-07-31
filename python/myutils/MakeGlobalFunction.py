import copy
def GlobalFunction(inputAll):
    (objectClass,functionName, inputFunction) = input_
    objectCopy = copy.deepcopy(objectClass)
    function = getattr(objectCopy,functionName)
    return function(*inputFunction)

