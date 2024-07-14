#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 09:38:58 2023

@author: carson16
"""

import numpy as np

import jax
import jax.numpy as jnp
jax.config.update("jax_enable_x64", True)

# Some more helper functions to go from voigt to matrix format and vice versa
def voigtNotation(mat):
    return jnp.asarray([mat[0, 0], mat[1, 1], mat[2, 2], mat[1, 2], mat[0, 2], mat[0, 1]])

def matNotation(voigt):
    return jnp.asarray([[voigt[0], voigt[5], voigt[4]],
                       [voigt[5], voigt[1], voigt[3]],
                       [voigt[4], voigt[3], voigt[2]]])

def matNotation_np(voigt):
    return np.asarray([[voigt[0], voigt[5], voigt[4]],
                       [voigt[5], voigt[1], voigt[3]],
                       [voigt[4], voigt[3], voigt[2]]])

def effectiveTerm(mat):
    term1 = mat[0, 0] - mat[1, 1]
    term2 = mat[1, 1] - mat[2, 2]
    term3 = mat[2, 2] - mat[0, 0]
    term4 = mat[1, 2] * mat[1, 2] \
          + mat[0, 2] * mat[0, 2] \
          + mat[0, 1] * mat[0, 1]

    return jnp.sqrt(0.5 * (term1 * term1 + term2 * term2 + term3 * term3 + 6.0 *  term4))

def symmetricMat(A):
    return 0.5 * (A + A.T)

def skewMat(A):
    return 0.5 * (A - A.T)

def skew_to_vec(W):
    return jnp.asarray([W[2, 1], W[0, 2], W[1, 0]])

def normalize(A):
    return jnp.linalg.norm(A) * A

def vec_dev_effective(A):
    norm = 0.0
    for ival in A:
        norm += ival * ival

    norm = jax.lax.cond(
        norm > 1e16,
        lambda: 1e16,
        lambda: norm
    )

    norm = jax.lax.cond(
        norm < 10.0 * jnp.finfo(jnp.float64).eps,
        lambda: 0.0,
        lambda: jnp.sqrt(norm)
    )

    return jnp.sqrt(2.0/3.0) * norm

def sym_mat_to_vec_dev(A):
    return jnp.asarray([ jnp.sqrt(0.5) * (A[0, 0] - A[1, 1]),
                         jnp.sqrt(1.0 / 6.0) * (2.0 * A[2, 2] - A[0, 0] - A[1, 1]),
                         jnp.sqrt(2.0) * A[1, 0],
                         jnp.sqrt(2.0) * A[2, 0],
                         jnp.sqrt(2.0) * A[2, 1]
                        ])

def sym_vec_to_vec_dev(svec):
    return jnp.asarray([ jnp.sqrt(0.5) * (svec[0] - svec[1]),
                         jnp.sqrt(1.0 / 6.0) * (2.0 * svec[2] - svec[0] - svec[1]),
                         jnp.sqrt(2.0) * svec[5],
                         jnp.sqrt(2.0) * svec[4],
                         jnp.sqrt(2.0) * svec[3]
                        ])

def sym_mat_to_sym_vec(A):
    vec_d = sym_mat_to_vec_dev(A)
    return jnp.c_[vec_d, jnp.sqrt(1.0 / 3.0) * jnp.trace(A)]

def inner_prod_sym_vec(vec1, vec2):
    return jnp.sum(vec1 * vec2) + 2.0 * jnp.sum(vec1[3:6] * vec2[3:6])

def dev_vec_to_sym_vec(dev_vec6):
    # dev_vec6 = [dev_vec_5, volume_term]
    t1 = jnp.sqrt(0.5) * dev_vec6[0]
    t2 = jnp.sqrt(1.0 / 6.0) * dev_vec6[1]
    
    return jnp.asarray([ t1 - t2,
                         -t1 - t2,
                         jnp.sqrt(2.0 / 3.0) * dev_vec6[1],
                         jnp.sqrt(0.5) * dev_vec6[4],
                         jnp.sqrt(0.5) * dev_vec6[3],
                         jnp.sqrt(0.5) * dev_vec6[2],
                         -jnp.sqrt(1.0 / 3.0) * dev_vec6[5]
                        ])

def sym_press_vec_to_sym_vec(svecm):
    return jnp.asarray([ svecm[0] - svecm[6],
                         svecm[1] - svecm[6],
                         svecm[2] - svecm[6],
                         svecm[3],
                         svecm[4],
                         svecm[5]
                        ])

def matrix_to_p_q(mat):
    P = sym_mat_to_vec_dev(symmetricMat(mat))
    Q = skew_to_vec(skewMat(mat))
    
    return (P, Q)

def axis_ang_to_quat(axis_map):
    sin_theta = jnp.sin(0.5 * axis_map[0])
    return jnp.asarray([ jnp.cos(axis_map[0] * 0.5),
                         sin_theta * axis_map[1],
                         sin_theta * axis_map[2],
                         sin_theta * axis_map[3]])

def quat_prod(quat1, quat2):
    q1 = quat1[0] * quat2[0] - quat1[1] * quat2[1] - quat1[2] * quat2[2] - quat1[3] * quat2[3]
    q2 = quat1[0] * quat2[1] + quat1[1] * quat2[0] + quat1[2] * quat2[3] - quat1[3] * quat2[2]
    q3 = quat1[0] * quat2[2] - quat1[1] * quat2[3] + quat1[2] * quat2[0] + quat1[3] * quat2[1]
    q4 = quat1[0] * quat2[3] + quat1[1] * quat2[2] - quat1[2] * quat2[1] + quat1[3] * quat2[0]
    
    return jnp.asarray([q1, q2, q3, q4])

def exp_map_to_quat(exp_map):
    # norm = 0.0
    # for ival in exp_map:
    #     norm += ival * ival

    # inorm = jax.lax.cond(
    #     norm < jnp.finfo(jnp.float64).eps,
    #     lambda: 0.0,
    #     lambda: 1.0 / jnp.sqrt(norm)
    # )
    
    # axis_ang = jax.lax.cond(
    #     norm < jnp.finfo(jnp.float64).eps,
    #     lambda: jnp.asarray([norm, 1.0, inorm * exp_map[1], inorm * exp_map[2]]),
    #     lambda: jnp.asarray([norm, inorm * exp_map[0], inorm * exp_map[1], inorm * exp_map[2]])
    # )
    
    # return axis_ang_to_quat(axis_ang)
    angle2 = jnp.dot(exp_map, exp_map)
    angle = jax.lax.cond(
        angle2 < jnp.finfo(jnp.float64).eps,
        lambda: jnp.finfo(jnp.float64).eps,
        lambda: jnp.sqrt(angle2)
    )
    small_scale = scale = 0.5 - angle2 / 48 + angle2 * angle2 / 3840
    large_scale = jnp.sin(angle / 2) / angle
    scale = jnp.where(angle <= 1e-3, small_scale, large_scale)
    return jnp.hstack([jnp.cos(angle / 2), scale * exp_map])

def quat_to_rmat(quat):
    x0sq = quat[0] * quat[0]

    x1sq = quat[1] * quat[1]
    x2sq = quat[2] * quat[2]
    x3sq = quat[3] * quat[3]

    x0x1 = quat[0] * quat[1]
    x0x2 = quat[0] * quat[2]
    x0x3 = quat[0] * quat[3]

    x1x2 = quat[1] * quat[2]
    x1x3 = quat[1] * quat[3]

    x2x3 = quat[2] * quat[3]
    
    return jnp.asarray([
        [x0sq + x1sq - x2sq - x3sq, 2.0 * (x1x2 - x0x3), 2.0 * (x1x3 + x0x2)],
        [2.0 * (x1x2 + x0x3), x0sq - x1sq + x2sq - x3sq, 2.0 * (x2x3 - x0x1)],
        [2.0 * (x1x3 - x0x2), 2.0 * (x2x3 + x0x1), x0sq - x1sq - x2sq + x3sq]
        ])

def update_quat_rot(dquat, quat_n):
    return quat_prod(quat_n, dquat)

def rot_mat_to_rot_mat5(rot_mat):
    rot_mat5 = jnp.zeros((5, 5)) 
    
    sqr3 = jnp.sqrt(3)
    
    rot_mat5 = rot_mat5.at[0, 0].set(0.5 * (rot_mat[0, 0] * rot_mat[0, 0] - rot_mat[0, 1] * rot_mat[0, 1] - rot_mat[1, 0] * rot_mat[1, 0] + rot_mat[1, 1] * rot_mat[1, 1]))
    rot_mat5 = rot_mat5.at[0, 1].set(sqr3 * 0.5 * (rot_mat[0, 2] * rot_mat[0, 2] - rot_mat[1, 2] * rot_mat[1, 2]))
    rot_mat5 = rot_mat5.at[0, 2].set(rot_mat[0, 0] * rot_mat[0, 1] - rot_mat[1, 0] * rot_mat[1, 1])
    rot_mat5 = rot_mat5.at[0, 3].set(rot_mat[0, 0] * rot_mat[0, 2] - rot_mat[1, 0] * rot_mat[1, 2])
    rot_mat5 = rot_mat5.at[0, 4].set(rot_mat[0, 1] * rot_mat[0, 2] - rot_mat[1, 1] * rot_mat[1, 2])
    rot_mat5 = rot_mat5.at[1, 0].set(sqr3 * 0.5 * (rot_mat[2, 0] * rot_mat[2, 0] - rot_mat[2, 1] * rot_mat[2, 1]))
    rot_mat5 = rot_mat5.at[1, 1].set(1.5 * rot_mat[2, 2] * rot_mat[2, 2] - 0.5)
    rot_mat5 = rot_mat5.at[1, 2].set(sqr3 * rot_mat[2, 0] * rot_mat[2, 1])
    rot_mat5 = rot_mat5.at[1, 3].set(sqr3 * rot_mat[2, 0] * rot_mat[2, 2])
    rot_mat5 = rot_mat5.at[1, 4].set(sqr3 * rot_mat[2, 1] * rot_mat[2, 2])
    rot_mat5 = rot_mat5.at[2, 0].set(rot_mat[0, 0] * rot_mat[1, 0] - rot_mat[0, 1] * rot_mat[1, 1])
    rot_mat5 = rot_mat5.at[2, 1].set(sqr3 * rot_mat[0, 2] * rot_mat[1, 2])
    rot_mat5 = rot_mat5.at[2, 2].set(rot_mat[0, 0] * rot_mat[1, 1] + rot_mat[0, 1] * rot_mat[1, 0])
    rot_mat5 = rot_mat5.at[2, 3].set(rot_mat[0, 0] * rot_mat[1, 2] + rot_mat[0, 2] * rot_mat[1, 0])
    rot_mat5 = rot_mat5.at[2, 4].set(rot_mat[0, 1] * rot_mat[1, 2] + rot_mat[0, 2] * rot_mat[1, 1])
    rot_mat5 = rot_mat5.at[3, 0].set(rot_mat[0, 0] * rot_mat[2, 0] - rot_mat[0, 1] * rot_mat[2, 1])
    rot_mat5 = rot_mat5.at[3, 1].set(sqr3 * rot_mat[0, 2] * rot_mat[2, 2])
    rot_mat5 = rot_mat5.at[3, 2].set(rot_mat[0, 0] * rot_mat[2, 1] + rot_mat[0, 1] * rot_mat[2, 0])
    rot_mat5 = rot_mat5.at[3, 3].set(rot_mat[0, 0] * rot_mat[2, 2] + rot_mat[0, 2] * rot_mat[2, 0])
    rot_mat5 = rot_mat5.at[3, 4].set(rot_mat[0, 1] * rot_mat[2, 2] + rot_mat[0, 2] * rot_mat[2, 1])
    rot_mat5 = rot_mat5.at[4, 0].set(rot_mat[1, 0] * rot_mat[2, 0] - rot_mat[1, 1] * rot_mat[2, 1])
    rot_mat5 = rot_mat5.at[4, 1].set(sqr3 * rot_mat[1, 2] * rot_mat[2, 2])
    rot_mat5 = rot_mat5.at[4, 2].set(rot_mat[1, 0] * rot_mat[2, 1] + rot_mat[1, 1] * rot_mat[2, 0])
    rot_mat5 = rot_mat5.at[4, 3].set(rot_mat[1, 0] * rot_mat[2, 2] + rot_mat[1, 2] * rot_mat[2, 0])
    rot_mat5 = rot_mat5.at[4, 4].set(rot_mat[1, 1] * rot_mat[2, 2] + rot_mat[1, 2] * rot_mat[2, 1])
    
    return rot_mat5

def mat35_da_A_oper_b_d(dev_a):
    # The opertor returned performs the following equivalent mapping
    # v = Ma
    # V = AB - BA
    # where a = deviatoric_vec(A), b = deviatoric_vec(B)
    # If we need the reverse operator aka BA - AB then
    # take negative of A
    # We can note that V is a skew matrix and therefore the opertor of
    # M transforms things from the 5d vec space to 3d vec space

    m35 = jnp.zeros((3, 5))

    m35 = m35.at[0, 0].set(dev_a[4] * 0.5)
    m35 = m35.at[1, 0].set(dev_a[3] * 0.5)
    m35 = m35.at[2, 0].set(-dev_a[2])

    m35 = m35.at[0, 1].set(dev_a[4] * 0.5 * jnp.sqrt(3))
    m35 = m35.at[1, 1].set(-dev_a[3] * 0.5 * jnp.sqrt(3))
    m35 = m35.at[2, 1].set(0.0)

    m35 = m35.at[0, 2].set(-dev_a[3] * 0.5)
    m35 = m35.at[1, 2].set(dev_a[4] * 0.5)
    m35 = m35.at[2, 2].set(dev_a[0])

    m35 = m35.at[0, 3].set(dev_a[2] * 0.5)
    m35 = m35.at[1, 3].set(0.5 * (jnp.sqrt(3) * dev_a[1] - dev_a[0]))
    m35 = m35.at[2, 3].set(-dev_a[4] * 0.5)

    m35 = m35.at[0, 4].set(-0.5 * (jnp.sqrt(3) * dev_a[1] + dev_a[0]))
    m35 = m35.at[1, 4].set(-dev_a[2] * 0.5)
    m35 = m35.at[2, 4].set(dev_a[3] * 0.5)
    
    return m35

def mtan_conv_sd_svec(mtanSD_vecds_raw, l_ddsdde_gamma):
    C = np.zeros((6,6))
    mtanSD = np.zeros((6,6))
    t1 = np.zeros(6)
    t2 = np.zeros(6)
    t3 = np.zeros(6)

    # mtanSD = T . mtanSD_vecds . T^{-1}
    # C = T . mtanSD_vecds
    # C(i,:) = T(i,k) . mtanSD_vecds(k,:) -- sum over k
    #
    t3 = mtanSD_vecds_raw[-1, :] / np.sqrt(3.0)
    t1 = mtanSD_vecds_raw[0, :] / np.sqrt(2.0)
    t2 = mtanSD_vecds_raw[1, :] / np.sqrt(6.0)

    C[0, :] = t1 - t2 + t3
    C[1, :] = -t1 - t2 + t3
    C[2, :] = np.sqrt(2.0 / 3.0) * mtanSD_vecds_raw[1, :] + t3
    C[3, :] = np.sqrt(1.0 / 2.0) * mtanSD_vecds_raw[4, :]
    C[4, :] = np.sqrt(1.0 / 2.0) * mtanSD_vecds_raw[3, :]
    C[5, :] = np.sqrt(1.0 / 2.0) * mtanSD_vecds_raw[2, :]

    # mtanSD = C . T^{-1}
    # mtanSD(:,j) = C(:,k) . [T^{-1}](k,j) -- sum over k
    #

    t3 = C[:, -1] / np.sqrt(3.0)
    t2 = C[:, 0] / np.sqrt(2.0)
    t1 = C[:, 1] / np.sqrt(6.0)

    mtanSD[:, 0] = t1 - t2 + t3
    mtanSD[:, 1] = -t1 - t2 + t3
    mtanSD[:, 2] = np.sqrt(2.0/3.0) * C[:, 1] + t3

    val = jax.lax.cond(
       l_ddsdde_gamma,
        lambda: 1.0 / np.sqrt(2.0),
        lambda: np.sqrt(2.0)
    )

    mtanSD[:, 3] = C[:, 4] * val
    mtanSD[:, 4] = C[:, 3] * val
    mtanSD[:, 5] = C[:, 2] * val

    return mtanSD
