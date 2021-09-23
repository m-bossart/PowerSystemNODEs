#indices you keep from full order parameter vector for reduced order paramter vector
#eg p_reduced = p_full[i__reduce_p], where p_full includes the transformer values and ss voltages
i__reduce_p = [1,2,3,4,5,6,7,8,9,10,12,13,16,24,25,26,27,28,29,30,21] 

#indices you keep from full order state vector for reduced order state vector
#eg u_reduced = u_full[i__reduce_u]
i__reduce_u = [3,6,8,11,12,13,14,16,17,20,21,22,23] 

#Indices of states in the surrogate for saving/ploting/calculating loss
i__ir_filter = 5    #what if this isn't a state? 
i__ii_filter = 19
i__ir_nn = 20
i__ii_nn = 21
i__ir_out = 22
i__ii_out = 23

