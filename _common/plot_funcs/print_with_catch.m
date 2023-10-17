function print_with_catch(filename,file_type)
    try
        print('-vector',file_type,filename); pause(0.2);
    catch msg
       disp(msg.message); 
    end
end