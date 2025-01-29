import numpy as np


def doubleAreaPoly(vertices):
    verticesLen = len(vertices)
    indices = np.roll(range(0, verticesLen), -1)
    return vertices[:, 0].T @ vertices[indices, 1] - \
        vertices[:, 1].T @ vertices[indices, 0]


def DPlaneStress(E, v):
    return E/(1-v ** 2) * np.array([[1, v, 0], [v, 1, 0], [0, 0, (1-v)/2]])


def BNatural2Global(BCoef, J):
    # B2
    invJ = np.linalg.inv(J)
    B2 = np.zeros((4, 4))
    B2[np.ix_([0, 1], [0, 1])] = invJ
    B2[np.ix_([2, 3], [2, 3])] = invJ

    # B3
    BCoefLen = len(BCoef[0, :])
    B3 = np.zeros((4, BCoefLen*2))
    B3[np.ix_([0, 1], range(0, BCoefLen*2, 2))] = BCoef
    B3[np.ix_([2, 3], range(1, BCoefLen*2, 2))] = BCoef

    return np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 1, 0]])@B2@B3


def BCoefLinearTri(vertices):
    b1 = vertices[1, 1]-vertices[2, 1]
    b2 = vertices[2, 1]-vertices[0, 1]
    b3 = vertices[0, 1]-vertices[1, 1]
    c1 = vertices[2, 0]-vertices[1, 0]
    c2 = vertices[0, 0]-vertices[2, 0]
    c3 = vertices[1, 0]-vertices[0, 0]
    return np.array([
        [b1, 0, b2, 0, b3, 0],
        [0, c1, 0, c2, 0, c3],
        [c1, b1, c2, b2, c3, b3]
    ])


def BCoefQuadraticTri(L1, L2):
    return np.array([
        [4*L1-1, 0, -3+4*L1+4*L2, 4*L2, -4*L2, 4-4*L2-8*L1],
        [0, 4*L2-1, -3+4*L1+4*L2, 4*L1, 4-4*L1-8*L2, -4*L1]
    ])


def stiffnessLinear(D, t, vertices, BCoefFunc):
    doubleArea = doubleAreaPoly(vertices)
    B = BCoefFunc(vertices) / doubleArea
    return B.T @ D @ B * t * doubleArea * 0.5, B


def stiffnessQuadratic(D, t, vertices, BCoefFunc):
    k = 0
    B = 0
    samplingPoints = np.array([[.5, .5], [.5, 0], [0, .5]])
    weights = [1/6, 1/6, 1/6]
    samplingPointsLen = len(samplingPoints[:, 0])
    for i in range(samplingPointsLen):
        BCoef = BCoefFunc(samplingPoints[i, 0], samplingPoints[i, 1])
        J = BCoef@vertices
        BGlob = BNatural2Global(BCoef, J)
        detJ = np.linalg.det(J)
        k += t*weights[i]*detJ*BGlob.T@D@BGlob
        B += BGlob

    B = B/samplingPointsLen
    return k, B


def evaluate(vertices, indices, forces, E, v, t, DoF, KConstructor, BCoefConstructor):
    indicesLen = len(indices)
    verticesLen = len(vertices)
    globalK = np.zeros((verticesLen * 2, verticesLen * 2))
    F = np.reshape(forces, (len(forces)*2, 1))
    D = DPlaneStress(E, v)
    e = np.zeros((indicesLen, 3))
    s = np.zeros(e.shape)
    Bs = np.zeros((3*indicesLen, 2*len(indices[0, :])))

    for i in range(indicesLen):
        verts = getFace(vertices, indices, i)
        K, B = KConstructor(D, t, verts, BCoefConstructor)

        # collect B for strain calculation
        # 0 1 2 3
        # 0 3 6 9 i * 3
        Bs[3*i:3*i+3, :] = B

        # assemble global stiffness matrix
        faceLen = len(indices[i, :])
        kIdx = np.zeros((faceLen*2, 1))
        kIdx[0:faceLen*2:2] = (2*indices[i, :]).reshape(faceLen, 1)
        kIdx[1:faceLen*2:2] = (2*indices[i, :]+1).reshape(faceLen, 1)
        kIdx = kIdx.astype(int)
        globalK[np.ix_(kIdx.flat, kIdx.flat)] += K

    # figure the row indices eliminated by DoF
    eliminateIdx = []
    dofLen = len(DoF)
    for i in range(dofLen):
        isHorizontal = i == 1
        for j in range(verticesLen):
            dof = DoF[i, j]
            if not dof:
                continue
            idx = 2 * j + isHorizontal * 1
            eliminateIdx.append(idx)

    # reduce matrices
    tempK = globalK.copy()
    tempIdx = range(len(globalK))
    tempK = np.delete(tempK, eliminateIdx, 0)
    tempK = np.delete(tempK, eliminateIdx, 1)
    tempIdx = np.delete(tempIdx, eliminateIdx)
    F = np.delete(F, eliminateIdx)

    # calculate displacement
    tempU = np.linalg.solve(tempK, F)
    u = np.zeros((2*verticesLen, 1))
    u[tempIdx] = tempU.reshape(len(tempU), 1)

    # calculate nodal forces
    F = globalK @ u

    # calculate element stresses and strains
    for i in range(indicesLen):
        faceLen = len(indices[i, :])
        uIdx = np.zeros((faceLen*2, 1))
        uIdx[0:faceLen*2:2] = (2*indices[i, :]).reshape(faceLen, 1)
        uIdx[1:faceLen*2:2] = (2*indices[i, :]+1).reshape(faceLen, 1)

        # element strain
        e[i, :] = (Bs[3*i:3*i+3, :] @ u[uIdx.astype(int).flat]
                   ).reshape(e[i, :].shape)

        # element stress
        s[i, :] = D@e[i, :].T

    # format output
    u = u.reshape(verticesLen, 2)
    F = F.reshape(verticesLen, 2)

    return u, F, e, s


def getFace(vertices, indices, faceIndex):
    return vertices[indices[faceIndex, :], :]
