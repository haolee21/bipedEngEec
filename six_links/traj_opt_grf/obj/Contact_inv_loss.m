function loss = Contact_inv_loss(x,p)

loss = 0;
loss = loss + Dyn_loss(x,p);
loss = loss + Gait_loss(x,p);
loss = loss + LCI_loss(x,p);



end