function plot_tracks(est,t1,t2)

  labelcount = countestlabels();
  colorarray = makecolorarray(labelcount);

  for k = t1:t2 % 1:length(est.X)
    if ~isempty(est.X{k})
        P= est.X{k}([1 3],:);
        for eidx = 1:size(P,2)
            yhline2 = line(P(1,eidx),P(2,eidx),'LineStyle','none','Marker','.','Markersize',15,'Color',colorarray.rgb(assigncolor(est.L{k}(:,eidx)),:));
        end
    end
  end
  
  
  function ca= makecolorarray(nlabels)
    lower= 0.1;
    upper= 0.9;
    rrr= rand(1,nlabels)*(upper-lower)+lower;
    ggg= rand(1,nlabels)*(upper-lower)+lower;
    bbb= rand(1,nlabels)*(upper-lower)+lower;
    ca.rgb= [rrr; ggg; bbb]';
    ca.lab= cell(nlabels,1);
    ca.cnt= 0;   
  end

  function idx= assigncolor(label)
    str= sprintf('%i*',label);
    tmp= strcmp(str,colorarray.lab);
    if any(tmp)
        idx= find(tmp);
    else
        colorarray.cnt= colorarray.cnt + 1;
        colorarray.lab{colorarray.cnt}= str;
        idx= colorarray.cnt;
    end
  end

  function count= countestlabels
    labelstack= [];
    for k=1:length(est.X)
        labelstack= [labelstack est.L{k}];
    end
    [c,~,~]= unique(labelstack','rows');
    count=size(c,1);
  end

end