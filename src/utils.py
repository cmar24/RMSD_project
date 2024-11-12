from scipy.spatial.transform import Rotation as R
import numpy as np
import matplotlib.pyplot as plt


def compute_rssd(true_atoms, p_atoms):
    """
    Compute the root sum square deviation between the true atoms and the predicted atoms
    :param true_atoms: atoms from the native structure
    :param p_atoms: atoms from the predictive structure
    :return: (float) root sum square deviation
    """
    rotation, rssd = R.align_vectors(true_atoms, p_atoms, return_sensitivity=False)
    return rotation, rssd


def plot_points(true_atoms, p_atoms, rotation):
    """
    Plot the atoms from the true structure, the predicted one, the rotated atoms that
    minimize the root sum square deviation and the final superimposed atoms
    :param true_atoms: atoms from the native structure
    :param p_atoms: atoms from the predictive structure
    :param rotation: output of R.align_vectors that contains the rotation matrix
    """
    rot_atoms = rotation.apply(p_atoms)
    # Compute translation
    translation = true_atoms.mean(axis=0) - rot_atoms.mean(axis=0)
    # Add to the rotated atoms
    trans_atoms = rot_atoms + translation
    # Plot the atoms
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(true_atoms[:, 0], true_atoms[:, 1], true_atoms[:, 2], label="True atoms")
    ax.scatter(p_atoms[:, 0], p_atoms[:, 1], p_atoms[:, 2], label="Predicted atoms")
    ax.scatter(rot_atoms[:, 0], rot_atoms[:, 1], rot_atoms[:, 2], label="Rotated atoms")
    ax.scatter(
        trans_atoms[:, 0],
        trans_atoms[:, 1],
        trans_atoms[:, 2],
        label="Translated atoms",
    )
    plt.legend()
    plt.show()


if __name__ == "__main__":
    points_A = np.array(
        [[4.267, -36.014, -5.602], [5.643, -34.581, -7.312], [7.753, -32.311, -6.707]]
    )
    points_B = np.array(
        [[-1.926, -18.243, -1.897], [-0.394, -18.282, -4.045], [0.418, -19.086, -7.040]]
    )
    rotation, rssd = compute_rssd(points_A, points_B)
    plot_points(points_A, points_B, rotation)
