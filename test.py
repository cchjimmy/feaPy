import fea as fea
import numpy as np
import matplotlib.pyplot as plt

# np.set_printoptions(suppress=True)
# help(np.set_printoptions)

E = 40
v = 0.3
t = 2
L = 40
F = np.array([0, -50])

# linear geometry
vertices = np.array([[0, 0], [L, 0], [60, 30], [L, 30], [0, 30]])
indices = np.array([[0, 3, 4], [0, 1, 3], [1, 2, 3]])

# global force vector
forces = np.zeros(vertices.shape)
forces[2, :] = F

# DoF, row order: horizontal, vertical
# 1 means DoF is eliminated. Indices correspond to vertices
DoF = np.zeros(vertices.shape).T
DoF[:, [0, 4]] = 1

u, F, e, s = fea.evaluate(vertices, indices, forces, E * 1000, v, t,
                          DoF, fea.stiffnessLinear, fea.BCoefLinearTri)

np.testing.assert_allclose(u, np.array(
    [[0, 0], [-.0013, -.0039], [.0012, -.0077], [.0013, -.0038], [0, 0]]), 0, 1e-4)

np.testing.assert_allclose(F, np.array(
    [[100, -10.5146], [0, 0], [0, -50], [0, 0], [-100, 60.5146]]), 0, 1e-4)

np.testing.assert_allclose(
    e, 1e-3*np.array([[0.0314, 0, -0.0952], [-.0327, .0044, -.0131], [-.0013, .0044, -.1083]]), 0, 1e-4)

np.testing.assert_allclose(s, np.array(
    [[1.3796, .4139, -1.4653], [-1.3796, -.2386, -.2013], [0, .1752, -1.6667]]), 0, 1e-4)


verticesQ = np.array([[0, 0], [L*.5, 0], [L, 0], [L+(60-L)*.5, 15], [60, 30],
                      [L+(60-L)*.5, 30], [L, 30], [L*.5, 30], [0, 30], [0, 15], [L*.5, 15], [L, 15]])
indicesQ = np.array(
    [[0, 6, 8, 10, 7, 9], [0, 2, 6, 1, 11, 10], [2, 4, 6, 3, 5, 11]])

DoF = np.zeros(verticesQ.shape).T
DoF[:, [0, 8]] = 1

forces = np.zeros(verticesQ.shape)
forces[4, :] = [0, -50]

uQ, FQ, eQ, sQ = fea.evaluate(verticesQ, indicesQ, forces, E*1000,
                              v, t, DoF, fea.stiffnessQuadratic, fea.BCoefQuadraticTri)

scale = 200
plt.subplot(1, 2, 1)
plt.scatter(vertices[:, 0], vertices[:, 1])
displaced = vertices + u * scale
plt.scatter(displaced[:, 0], displaced[:, 1])

plt.subplot(1, 2, 2)
plt.scatter(verticesQ[:, 0], verticesQ[:, 1])
displaced = verticesQ+uQ*scale
plt.scatter(displaced[:, 0], displaced[:, 1])
plt.show()
