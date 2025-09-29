import smtplib
from email.mime.text import MIMEText
# 1. 导入 formataddr
from email.utils import formataddr
from email.header import Header
import argparse

def send_email(status, recipient_email):
    # 邮件服务器配置
    smtp_server = 'smtp.qq.com'
    smtp_port = 465
    sender_email = '1425623506@qq.com'
    sender_password = 'lmhisntzvkkagegi' # 你的邮箱授权码

    # 邮件内容
    if status == 'start':
        subject = '程序运行开始通知'
        body = '您好，程序已开始运行。'
    elif status == 'end':
        subject = '程序运行结束通知'
        body = '您好，程序已运行结束。'
    else:
        raise ValueError("Invalid status. Use 'start' or 'end'.")

    # 构造邮件
    message = MIMEText(body, 'plain', 'utf-8')
    
    # 2. 使用 formataddr 来格式化 From 和 To 字段
    # 这样可以确保格式绝对正确
    message['From'] = formataddr(('通知系统', sender_email))
    message['To'] = formataddr(('用户', recipient_email))
    
    # Subject 字段只包含文本，继续使用 Header 即可
    message['Subject'] = Header(subject, 'utf-8')

    server = None
    try:
        # 使用 smtplib.SMTP_SSL() 创建SSL加密连接
        server = smtplib.SMTP_SSL(smtp_server, smtp_port, timeout=30)
        server.login(sender_email, sender_password)
        server.sendmail(sender_email, [recipient_email], message.as_string())
        print(f"邮件已成功发送至 {recipient_email}")
        
    except smtplib.SMTPException as e:
        print(f"发送邮件失败，SMTP错误: {e}")
    except Exception as e:
        print(f"发送邮件失败，发生未知错误: {e}")
    finally:
        if server:
            server.quit()

# 示例调用 (这部分保持不变)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('status', choices=['start', 'end'], help="程序状态: 'start' 或 'end'")
    args = parser.parse_args()
    status = args.status
    recipient_email = '2010307209@stu.pku.edu.cn'
    send_email(status, recipient_email)