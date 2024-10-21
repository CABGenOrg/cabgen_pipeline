from subprocess import run
from os import getenv, path
from dotenv import load_dotenv

load_dotenv()

sender_email = getenv("SENDER_EMAIL") or ""
template_email_path = getenv("TEMPLATE_EMAIL_PATH") or ""


def send_email(recipient_email: str, subject: str, template: str):
    try:
        template_email = path.join(template_email_path, template)
        command_line = (f"mailx -r {sender_email} -s '{subject}' "
                        f"{recipient_email} < {template_email}")
        run(command_line, shell=True)
    except Exception as e:
        raise Exception(f"Failed to send e-mail.\n\n{e}")
