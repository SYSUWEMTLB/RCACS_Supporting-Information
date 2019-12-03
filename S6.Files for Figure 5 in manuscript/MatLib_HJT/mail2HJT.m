function mail2HJT(subject,content,varargin) 
if nargin==3
    file = varargin{1};
end
MailAddress = 'diwangeyong@163.com'; 
password = 'q8wyzl3.1415926';   
setpref('Internet','E_mail',MailAddress); 
setpref('Internet','SMTP_Server','smtp.163.com'); 
setpref('Internet','SMTP_Username',MailAddress); 
setpref('Internet','SMTP_Password',password); 
props = java.lang.System.getProperties; 
props.setProperty('mail.smtp.auth','true'); 
if nargin==3
    sendmail('jiatanghu@126.com ',subject,content,file);
else
    sendmail('jiatanghu@126.com ',subject,content);
end