import QTCIHaldane

δm = 1e-4
R = 10
tolerance = 1e-4
datadirectory = "."

evalooserror = false

QTCIHaldane.evaluatechern_haldane(δm, R; tolerance, datadirectory, evalooserror)
