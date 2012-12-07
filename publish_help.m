% Just a simple script to publish all the help files.
    delete('html/*');
    d = dir('*.m');
    for k = 1 : length(d)
        f = d(k).name;

        if strfind(f, 'example')
            publish(f);
        else
            publish(f, 'evalCode', false);
        end
    end
        

