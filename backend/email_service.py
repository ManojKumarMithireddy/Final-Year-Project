"""
email_service.py
─────────────────
Sends transactional emails via the Brevo (formerly Sendinblue) API.

Required environment variables (add to backend/.env):
    BREVO_API_KEY       → Your Brevo API key (Settings → API Keys)
    BREVO_SENDER_EMAIL  → Verified sender email address
    BREVO_SENDER_NAME   → (optional) Display name, defaults to "BioQuantum"
"""

import os
import sib_api_v3_sdk
from sib_api_v3_sdk.rest import ApiException


def _get_api_instance():
    configuration = sib_api_v3_sdk.Configuration()
    configuration.api_key["api-key"] = os.getenv("BREVO_API_KEY", "")
    return sib_api_v3_sdk.TransactionalEmailsApi(
        sib_api_v3_sdk.ApiClient(configuration)
    )


def send_verification_email(to_email: str, to_name: str, verify_url: str) -> None:
    """
    Send an HTML verification email to the newly registered user.

    Args:
        to_email:   Recipient email address.
        to_name:    Recipient display name.
        verify_url: The full verification URL containing the signed token.

    Raises:
        RuntimeError: If BREVO_API_KEY or BREVO_SENDER_EMAIL are not configured,
                      or if the Brevo API call itself fails.
    """
    api_key = os.getenv("BREVO_API_KEY", "")
    sender_email = os.getenv("BREVO_SENDER_EMAIL", "")
    sender_name = os.getenv("BREVO_SENDER_NAME", "BioQuantum")

    if not api_key or not sender_email:
        raise RuntimeError(
            "Brevo is not configured. Set BREVO_API_KEY and BREVO_SENDER_EMAIL in backend/.env"
        )

    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="UTF-8" />
      <meta name="viewport" content="width=device-width, initial-scale=1.0" />
      <title>Verify your BioQuantum account</title>
    </head>
    <body style="margin:0;padding:0;background:#0f172a;font-family:Inter,Arial,sans-serif;">
      <table width="100%" cellpadding="0" cellspacing="0" style="padding:40px 0;">
        <tr>
          <td align="center">
            <table width="520" cellpadding="0" cellspacing="0"
                   style="background:#1e293b;border-radius:16px;padding:40px;border:1px solid #334155;">
              <!-- Logo / wordmark -->
              <tr>
                <td align="center" style="padding-bottom:24px;">
                  <div style="display:inline-block;background:linear-gradient(135deg,#3b82f6,#7c3aed);
                              border-radius:12px;padding:12px 20px;">
                    <span style="color:#fff;font-size:22px;font-weight:700;letter-spacing:0.5px;">
                      🧬 BioQuantum
                    </span>
                  </div>
                </td>
              </tr>
              <!-- Heading -->
              <tr>
                <td align="center" style="padding-bottom:12px;">
                  <h1 style="color:#f1f5f9;font-size:22px;margin:0;font-weight:600;">
                    Confirm your email address
                  </h1>
                </td>
              </tr>
              <!-- Body text -->
              <tr>
                <td align="center" style="padding-bottom:28px;">
                  <p style="color:#94a3b8;font-size:15px;line-height:1.6;margin:0;">
                    Hi {to_name}, welcome to BioQuantum!<br/>
                    Click the button below to verify your email and activate your account.<br/>
                    This link expires in <strong style="color:#cbd5e1;">1 hour</strong>.
                  </p>
                </td>
              </tr>
              <!-- CTA button -->
              <tr>
                <td align="center" style="padding-bottom:28px;">
                  <a href="{verify_url}"
                     style="display:inline-block;background:linear-gradient(135deg,#3b82f6,#7c3aed);
                            color:#fff;font-size:15px;font-weight:600;text-decoration:none;
                            padding:14px 36px;border-radius:10px;letter-spacing:0.3px;">
                    Verify my email →
                  </a>
                </td>
              </tr>
              <!-- Fallback link -->
              <tr>
                <td align="center" style="padding-bottom:16px;">
                  <p style="color:#64748b;font-size:12px;margin:0;">
                    Or copy and paste this link:<br/>
                    <a href="{verify_url}" style="color:#7c3aed;word-break:break-all;">
                      {verify_url}
                    </a>
                  </p>
                </td>
              </tr>
              <!-- Footer note -->
              <tr>
                <td align="center">
                  <p style="color:#475569;font-size:11px;margin:0;">
                    If you didn't create a BioQuantum account, you can safely ignore this email.
                  </p>
                </td>
              </tr>
            </table>
          </td>
        </tr>
      </table>
    </body>
    </html>
    """

    send_smtp_email = sib_api_v3_sdk.SendSmtpEmail(
        to=[{"email": to_email, "name": to_name}],
        sender={"email": sender_email, "name": sender_name},
        subject="Verify your BioQuantum email address",
        html_content=html_content,
    )

    try:
        api_instance = _get_api_instance()
        api_instance.send_transac_email(send_smtp_email)
    except ApiException as e:
        raise RuntimeError(f"Brevo API error: {e}") from e


def send_resend_verification_email(to_email: str, to_name: str, verify_url: str) -> None:
    """Convenience wrapper — identical email but with a slightly different subject line."""
    send_verification_email(to_email, to_name, verify_url)


def send_password_reset_email(to_email: str, to_name: str, reset_url: str) -> None:
    """
    Send a password-reset email with a signed reset link.

    Args:
        to_email:  Recipient email address.
        to_name:   Recipient display name.
        reset_url: The full reset URL containing the signed token.

    Raises:
        RuntimeError: If Brevo credentials are missing or the API call fails.
    """
    api_key = os.getenv("BREVO_API_KEY", "")
    sender_email = os.getenv("BREVO_SENDER_EMAIL", "")
    sender_name = os.getenv("BREVO_SENDER_NAME", "BioQuantum")

    if not api_key or not sender_email:
        raise RuntimeError(
            "Brevo is not configured. Set BREVO_API_KEY and BREVO_SENDER_EMAIL in backend/.env"
        )

    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="UTF-8" />
      <meta name="viewport" content="width=device-width, initial-scale=1.0" />
      <title>Reset your BioQuantum password</title>
    </head>
    <body style="margin:0;padding:0;background:#0f172a;font-family:Inter,Arial,sans-serif;">
      <table width="100%" cellpadding="0" cellspacing="0" style="padding:40px 0;">
        <tr>
          <td align="center">
            <table width="520" cellpadding="0" cellspacing="0"
                   style="background:#1e293b;border-radius:16px;padding:40px;border:1px solid #334155;">
              <!-- Logo / wordmark -->
              <tr>
                <td align="center" style="padding-bottom:24px;">
                  <div style="display:inline-block;background:linear-gradient(135deg,#3b82f6,#7c3aed);
                              border-radius:12px;padding:12px 20px;">
                    <span style="color:#fff;font-size:22px;font-weight:700;letter-spacing:0.5px;">
                      🧬 BioQuantum
                    </span>
                  </div>
                </td>
              </tr>
              <!-- Heading -->
              <tr>
                <td align="center" style="padding-bottom:12px;">
                  <h1 style="color:#f1f5f9;font-size:22px;margin:0;font-weight:600;">
                    Reset your password
                  </h1>
                </td>
              </tr>
              <!-- Body text -->
              <tr>
                <td align="center" style="padding-bottom:28px;">
                  <p style="color:#94a3b8;font-size:15px;line-height:1.6;margin:0;">
                    Hi {to_name},<br/>
                    We received a request to reset your BioQuantum password.<br/>
                    Click the button below to choose a new password.<br/>
                    This link expires in <strong style="color:#cbd5e1;">1 hour</strong>.
                  </p>
                </td>
              </tr>
              <!-- CTA button -->
              <tr>
                <td align="center" style="padding-bottom:28px;">
                  <a href="{reset_url}"
                     style="display:inline-block;background:linear-gradient(135deg,#3b82f6,#7c3aed);
                            color:#fff;font-size:15px;font-weight:600;text-decoration:none;
                            padding:14px 36px;border-radius:10px;letter-spacing:0.3px;">
                    Reset my password →
                  </a>
                </td>
              </tr>
              <!-- Fallback link -->
              <tr>
                <td align="center" style="padding-bottom:16px;">
                  <p style="color:#64748b;font-size:12px;margin:0;">
                    Or copy and paste this link:<br/>
                    <a href="{reset_url}" style="color:#7c3aed;word-break:break-all;">
                      {reset_url}
                    </a>
                  </p>
                </td>
              </tr>
              <!-- Footer note -->
              <tr>
                <td align="center">
                  <p style="color:#475569;font-size:11px;margin:0;">
                    If you didn't request a password reset, you can safely ignore this email.
                    Your password will not be changed.
                  </p>
                </td>
              </tr>
            </table>
          </td>
        </tr>
      </table>
    </body>
    </html>
    """

    send_smtp_email = sib_api_v3_sdk.SendSmtpEmail(
        to=[{"email": to_email, "name": to_name}],
        sender={"email": sender_email, "name": sender_name},
        subject="Reset your BioQuantum password",
        html_content=html_content,
    )

    try:
        api_instance = _get_api_instance()
        api_instance.send_transac_email(send_smtp_email)
    except ApiException as e:
        raise RuntimeError(f"Brevo API error: {e}") from e
