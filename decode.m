function msg_rx = decode(data_in)
% Tim Drysdale
% 19 Nov 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%code properties (fixed in the encode function)
n0 = 2;
k0 = 1;
L =  3;
G = [ 1, 0, 1;
      1, 1, 1;];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up variables with hardcoded trellis properties (worked out by hand)

%state assignment
a = 1;
b = 2;
c = 3;
d = 4;

output =   [ 0,0; ...  %in branch order
             1,1; ...
             1,1; ...
             0,0; ...
             0,1; ...
             1,0; ...
             1,0; ...
             0,1; ];
          
from =  [   a; ...   %from state, in branch order
            a; ...
            b; ...
            b; ...
            c; ...
            c; ...
            d; ...
            d; ];
        
to  =  [    a; ...   %to state, in branch order
            c; ...
            a; ...
            c; ...
            b; ...
            d; ...
            b; ...
            d; ]; 

td = round(length(data_in))./n0; %trellis depth

sequence = reshape(data_in,n0,[])';

B = zeros(8,td);  %branch metric, branch order

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the branch metrics (this depends only on output and sequence)
for level = 1:td
    for ii=1:8
        B(ii,level) = hamming(output(ii,:),sequence(level,:));
    end
end

G = zeros(8,td);   %path metric, branch order
P = zeros(8,td+1); %path metric, state order
S = zeros(8,td+1); %survivor, branch order

P(:,1) = [      0; ...  %must start at state a, by definition of code
                0; ...
                99; ...
                99; ...
                99; ...
                99; ...
                99; ...
                99;];
            
S(:,1) = [      0; ... %no survivors before 1st state, by definition
                0; ...
                0; ...
                0; ...
                0; ...
                0; ...
                0; ...
                0;];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate path metrics and survivor branchs
for level = 1:td
    G(:,level)  = B(:,level) + P(:,level);          % branch order
    P(:,level+1)= [ min(G(1,level),G(3,level)); ... % state order
                    min(G(1,level),G(3,level)); ...
                    min(G(5,level),G(7,level)); ... 
                    min(G(5,level),G(7,level)); ...
                    min(G(2,level),G(4,level)); ... 
                    min(G(2,level),G(4,level)); ...
                    min(G(6,level),G(8,level)); ... 
                    min(G(6,level),G(8,level)); ]; 
    
    S(:,level+1)= [ G(1,level)<=G(3,level); ...  %branch order
                    G(2,level)<=G(4,level); ... 
                    G(3,level)<=G(1,level); ...
                    G(4,level)<=G(2,level); ...
                    G(5,level)<=G(7,level); ...  
                    G(6,level)<=G(8,level); ... 
                    G(7,level)<=G(5,level); ...
                    G(8,level)<=G(6,level);];

end % for 1:td

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%trace back - note that trace back does not check G, but it might need to??

W = zeros(1,td+1);

%find path to final state
if xor(S(1,td+1),S(3,td+1)) == 0
     %error('Two surviving paths')
     W(td+1) = 1;
end

if S(1,td+1) == 1
    W(td+1) = 1;
else
    W(td+1) = 3;
end

%now do rest of the levels (L)
for L = td:-1:1
    switch from(W(L+1))
          case a, 
            if S(1,L) == 1
                W(L) = 1;
            else
                W(L) = 3;
            end
          case b,
            if S(5,L) == 1
                W(L) = 5;
            else
                W(L) = 7;
            end
        case c,
            if S(2,L) == 1
                W(L) = 2;
            else
                W(L) = 4;
            end
        case d,
            if S(6,L) == 1
                W(L) = 6;
            else
                W(L) = 8;
            end
    end % switch
end %for L
    
state_out = [1;to(W)];

data_from_state = [ 0; ...
                    0; ...
                    1; ...
                    1; ];

decoded_data = data_from_state(state_out);

msg_rx = decoded_data(3:end-2)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %visualise the trellis
% %td = td-2;
% Tedge = zeros(4,4,td+1);
% 
% Tedge(a,a,2) = 1;
% Tedge(a,c,2) = 1;
% 
% Tedge(a,a,3) = 1;
% Tedge(a,c,3) = 1;
% Tedge(c,b,3) = 1;
% Tedge(c,d,3) = 1;
% 
% Tedge(a,a,4:td-1) = 1;
% Tedge(a,c,4:td-1) = 1;
% Tedge(b,a,4:td-1) = 1;
% Tedge(b,c,4:td-1) = 1;
% Tedge(c,b,4:td-1) = 1;
% Tedge(c,d,4:td-1) = 1;
% Tedge(d,b,4:td-1) = 1;
% Tedge(d,d,4:td-1) = 1;
% 
% Tedge(a,a,4:td-0) = 1;
% Tedge(b,a,4:td-0) = 1;
% Tedge(c,b,4:td-0) = 1;
% Tedge(d,b,4:td-0) = 1;
% 
% Tedge(a,a,td+1) = 1;
% Tedge(b,a,td+1) = 1;
% 
% %figure()
% for ii=1:td+1
%     [y1,y2] = find(Tedge(:,:,ii)==1);
%     for jj=1:length(y1)
%         h = line([ii-1,ii],[4-y1(jj),4-y2(jj)]);
%         set(h,'Color','black','LineWidth',[0.5],'LineStyle',':');
%     end
% end
% 
% %[X,Y] = meshgrid([1:td+1],[0:3]);
% %hold on;
% %plot(X,Y,'k.','Markersize',20);
% 
% %plot([1:td+1],4-state_out(2:end),'rs','Markersize',20)
% 
% for ii=2:td+1
%     y1 = 4-state_out(ii);
%     y2 = 4-state_out(ii+1);
%     h = line([ii-1,ii],[y1,y2]);
%     set(h,'Color','red','LineWidth',[1.0])
% end
% 
% %tidy up the labels
% %ax = gca;
% %set(ax,'YTickLabel',{'d/11','c/10','b/01','a/00'})
% %set(ax,'YTick',[0,1,2,3])
% %ylabel('State')
% %xlabel('time')
% %set(ax,'XTick',[0:td+1])