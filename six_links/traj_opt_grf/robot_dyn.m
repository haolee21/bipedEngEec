function dydt=robot_dyn(t,x,text,u,f_toe,f_heel,p)
u = interp1(text,u.',t);
f_toe = interp1(text,f_toe.',t);
f_heel = interp1(text,f_heel.',t);


f = f_dyn(x,p,f_toe.',f_heel.',u.');


dydt = f;

end