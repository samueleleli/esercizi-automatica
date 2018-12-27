%%by jacopo ferretti
%%17/12/18
%%enjoy!
clear all; clc;
format rat;

%%% INSERISCI QUI LE TUE MATRICI %%%
A=[-3 -1/2 -2;-10 1/2 -5;4 1 3];
B=[-1;-2;2];
C=[1 0 1];


%%%come avete notato, valgono solo le 3x3%%
fprintf('se non lo hai ancora fatto, inserisci nell editor:')
fprintf('\nla matrice A nel formato [1 2 3;4 5 6;7 8 9]');
fprintf('\nla matrice B nel formato [1;2;3]');
fprintf('\nla matrice C nel formato [1 2 3]');
fprintf('\n\nper salvare i dati premi su run');
fprintf('\n\npress any key\n');

pause;

%l'utente ha inserito tuttle le matrici
%calcolamo la matrice R di raggiungibilità

fprintf('\nla matrice della raggiungibilità è R:\n')
R=[B A*B A*A*B];
disp(R);
%calcolamo la matrice O di osservabilità

fprintf('\nla matrice dell osservabilità è O:\n')
O=[C;C*A;C*A*A];
disp(O);

%%trovo Xr

Xr=R(:,1);  %%inizializzo Xr con la prima colonna

%%inizializzo i e j, sono dei segnalini per controllare se i vettori sono
%%linearmente indiendenti
i=1;
j=2;
k=0;

while rank(Xr)<rank(R) && k~=3
    
    if size(Xr,2)~=rank([Xr R(:,j)])
        Xr=[Xr R(:,j)];
           
    end
    %%questo è per cambiare i le colonne delle matrici        
    if j==3
        i=2;
    end   
    
    if i==1
        j=3;
    end
    %%fine coso per cambiare le colonne
    
    k=k+1; %%controllo contro il loop infinito
end

if rank(Xr)==3             %%se il rango di Xr è 3 allora posso vederlo
    Xr=[1 0 0;0 1 0;0 0 1]; %%come una base canonica 
end

fprintf('\nla matrice Xr è\n');
disp(Xr);


%%trovaimo la matrice Xno

Xno=null(O,'r');

fprintf('\nla matrice Xno è\n');
disp(Xno);

%%trovo la matrice T
%%trovo la intersezione tra Xr e Xno, cioè T1

preintersezione=null([Xr -Xno],'r'); %%preliminare per fare l'intersezione
T1=[];                                    %%così trovo alfa e beta..                     


for j=1:size(preintersezione,2)

    alfa=preintersezione(1:size(Xr,2),j);
    
    T1=[T1 Xr*alfa];
end
T=T1;

% trovo t2 che competa ad una base di Xr
T2=[];
for j=1:size(Xr,2)
    
    if isempty(null([T Xr(:,j)]))
    T=[T Xr(:,j)];
    T2=[T2 Xr(:,j)];
    end
end

%%trovo T3 che completa ad una base di Xno                  

T3=[];
for j=1:size(Xno,2)
    
    if isempty(null([T Xno(:,j)]))
    T=[T Xno(:,j)];
    T3=[T3 Xno(:,j)];
    end
end

%%trovo T4 che completa ad una base di R3
I=eye(3);                  
T4=[];
for j=1:3
    
    if isempty(null([T I(:,j)]))
    T=[T I(:,j)];
    T4=[T4 I(:,j)];
    end
end 
%%% ho trovato la matrice T=[T1 T2 T3 T4]
fprintf('\nla matrice T è:\n')
disp(T)

%trovo A B C segneti

fprintf('\nla matrice A segnata è:\n')
Asegnato=inv(T)*A*T;
disp(Asegnato)
         
fprintf('\nla matrice B segnata è:\n')
Bsegnato=inv(T)*B;
disp(Bsegnato)

fprintf('\nla matrice c segnata è:\n')
Csegnato=C*T;
disp(Csegnato)
                                   
%%trovo gli autovaliri e gli autovettori                                    

autov=eig(Asegnato);

k=1; h=0; A11=[]; A22=[]; A33=[]; A44=[];
q1=size(T1,2); q2=size(T2,2); 
q4=size(T4,2); q3=size(T3,2);

%%colloco ogni autovalore al proprio sottosistema
%%scrivo gli autovalori in ordine e come sono

fprintf('\ngli autovalori sono :\n\n\n')
if ~isempty(T1)
    [A11,B1,C1,k]=calcolo_matrici(Asegnato,Bsegnato,Csegnato,k,q1);
    fprintf('raggiungibili non osservabili\n\n')
    h=calcolo_autovalori(h,A11,autov);
else
    clear T1;
end

if ~isempty(T2)
    [A22,B2,C2,k]=calcolo_matrici(Asegnato,Bsegnato,Csegnato,k,q2);    
    fprintf('\n\nraggiungibili osservabili\n\n')
    h=calcolo_autovalori(h,A22,autov);
else
    clear T2;
end

if ~isempty(T3)
    [A33,B3,C3,k]=calcolo_matrici(Asegnato,Bsegnato,Csegnato,k,q3);
    fprintf('\n\nnon raggiungibili non osservabili\n\n')
    h=calcolo_autovalori(h,A11,autov);
else
    clear T3;
end

if ~isempty(T4)
    [A44,B4,C4,k]=calcolo_matrici(Asegnato,Bsegnato,Csegnato,k,q4);
    fprintf('\n\nnon raggiungibili osservabili\n\n')
    h=calcolo_autovalori(h,A44,autov);
else
    clear T4;
end

%%trovo gli autovettori

for j=1:3
    
   fprintf('\n\nl autovettore di A associato a Lambda%d=%s è:\n',...
   j,strtrim(rats(autov(j,1))))                             %%sono andato accapo

   disp(T*(null((autov(j,1)*I)-Asegnato,'r')))
end

%%%trovo G(s) per la risposta
if ~isempty(A22)
    s=sym('s');
    Gs=C2*1/((s*eye(size(A22,2)))-A22)*B2;
    fprintf('G(s) per il calcolo della risposta è:\n\n')
    disp(Gs)
else
    fprintf(['non posso trovare G(s) perchè il sistema raggiungibile',...
    'e raggiungibile è vuoto'])
end

%%faccio un po' di pulizia per rendere più leggibili i risulati nella workspace

clear q1 q2 q3 q4 i j k s alfa I preintersezione sizepreint autov h;

function [Ax,By,Cz,k]=calcolo_matrici(Asegnato,Bsegnato,Csegnato,k,q)
    Ax=Asegnato(k:k-1+q,k:k-1+q);
    By=Bsegnato(k:k-1+q);
    Cz=Csegnato(k:k-1+q);
    k=k+q;
end

function k=calcolo_autovalori(k,I,autov)
    for j=1:size(I,2)
       if autov(j+k)<=0
           fprintf('%s convergente\n',strtrim(rats(autov(j+k))))%%fix per mandare 
        else                                                 %%in output il numero 
           fprintf('%s divergente\n',strtrim(rats(autov(j+k))))%%in forma razionale
        end
    end
    k=j;
end
