function fob=funobj(Pr, Eb0, T_EDF, Pc, Pp, Cb, Prmax, Prmin, Te)
            lamda=1e6;
            if(isrow(Pr))
                Pr=Pr';
            end
            Pb=Pr+Pp-Pc;
            [~,Pbr]=systeme(Pb, Eb0,Cb,Te);
            Prr=Pbr+Pc-Pp;
            Cm=lamda*((Prr<Prmin).*(Prmin-Prr)+(Prr>Prmax).*(Prr-Prmax));
            fob=sum(Cm+(T_EDF.*Prr)*Te/60);
end