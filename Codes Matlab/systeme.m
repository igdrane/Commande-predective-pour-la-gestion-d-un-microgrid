function [Eb, Pbr]=systeme(Pb, Eb0, Cb, Te)
            Ebmin=0.1*Cb;
            Ebmax=Cb;
            Pbr=Pb;
            N=length(Pb);
            Eb=zeros(N+1,1);
            Eb(1)=Eb0;
            for k=1:N
                Eb(k+1)=Eb(k)+Pb(k)*Te/60;
                if(Eb(k+1)<Ebmin)
                    Eb(k+1)=Ebmin;
                    Pbr(k)=(Ebmin-Eb(k))*60/Te;
                end
                if(Eb(k+1)>Ebmax)
                    Eb(k+1)=Ebmax;
                    Pbr(k)=(Ebmax-Eb(k))*60/Te;
                end
            end
            Eb=Eb(2:N+1);         
end