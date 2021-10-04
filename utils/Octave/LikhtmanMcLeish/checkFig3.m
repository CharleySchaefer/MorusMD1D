% Reproduces Fig 3 in Likhtman-McLeish, macromolecules, 2002
function checkFig3()
  Zrow=[1:1:1000];
  for i=1:length(Zrow)
    Z=Zrow(i);
    T(i)=get_renormalisation_reptation_time(Z);
    G(i)=get_renormalisation_elastic_modulus(Z);

  end

  figure
  plot(1./sqrt(Zrow),T); hold on;
  plot(1./sqrt(Zrow),G)
end 
