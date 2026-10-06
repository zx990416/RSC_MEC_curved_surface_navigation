function ang2 = ang_mod_360(ang1)
ang1 = ang1;
ang1(ang1<0) = ang1(ang1<0)+360;
ang1(ang1<0) = ang1(ang1<0)+360;
ang1(ang1<0) = ang1(ang1<0)+360;
ang1(ang1>360) = ang1(ang1>360)-360;
ang1(ang1>360) = ang1(ang1>360)-360;
ang1(ang1>360) = ang1(ang1>360)-360;
ang2 = ang1;
end
