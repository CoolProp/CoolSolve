# Deploying CoolSolve GUI on Ubuntu 24.04 with Apache Reverse Proxy

This guide covers deploying CoolSolve GUI on an Ubuntu 24.04 server behind Apache as a reverse proxy, with user authentication, automatic startup, and proper security.

## Architecture

```
Internet → Apache (port 443/80) → reverse proxy → CoolSolve (localhost:8550)
```

Apache handles TLS, authentication, and public access. CoolSolve listens only on `127.0.0.1:8550` and is not directly exposed to the internet.

---

## 1. Prerequisites

```bash
# Update system
sudo apt update && sudo apt upgrade -y

# Install Apache and required modules
sudo apt install -y apache2 apache2-utils

# Enable required Apache modules
sudo a2enmod proxy proxy_http proxy_wstunnel headers rewrite ssl auth_basic
sudo systemctl restart apache2
```

CoolSolve is assumed to be already compiled at `/home/coolsolve/CoolSolve/build/coolsolve`.

---

## 2. Dedicated System User

If not already done, create a dedicated system user for CoolSolve (no login shell, no home directory login):

```bash
# If the 'coolsolve' user doesn't exist yet:
sudo useradd -r -s /usr/sbin/nologin -d /home/coolsolve -m coolsolve

# Ensure ownership of the CoolSolve installation
sudo chown -R coolsolve:coolsolve /home/coolsolve/CoolSolve
```

If the user already exists with a normal shell, you can restrict it:

```bash
sudo usermod -s /usr/sbin/nologin coolsolve
```

---

## 3. File Permissions

```bash
# Set restrictive permissions on the installation directory
sudo chmod 750 /home/coolsolve/CoolSolve
sudo chmod 750 /home/coolsolve/CoolSolve/build

# The coolsolve binary must be executable
sudo chmod 755 /home/coolsolve/CoolSolve/build/coolsolve

# Examples directory must be readable by the coolsolve user
sudo chmod -R 755 /home/coolsolve/CoolSolve/examples

# Create a writable temp directory for session data
sudo mkdir -p /tmp/coolsolve_sessions
sudo chown coolsolve:coolsolve /tmp/coolsolve_sessions
sudo chmod 700 /tmp/coolsolve_sessions

# Recreate after reboot (/tmp is wiped; the service needs this path to exist)
echo 'd /tmp/coolsolve_sessions 0700 coolsolve coolsolve -' | sudo tee /etc/tmpfiles.d/coolsolve.conf
sudo systemd-tmpfiles --create /etc/tmpfiles.d/coolsolve.conf
```

---

## 4. Systemd Service (Auto-start)

Create the service file:

```bash
sudo tee /etc/systemd/system/coolsolve.service > /dev/null << 'EOF'
[Unit]
Description=CoolSolve GUI Server
After=network.target

[Service]
Type=simple
User=coolsolve
Group=coolsolve
WorkingDirectory=/home/coolsolve/CoolSolve/build
ExecStart=/home/coolsolve/CoolSolve/build/coolsolve --gui 8550 --no-browser
Restart=on-failure
RestartSec=5

# Security hardening
NoNewPrivileges=true
ProtectSystem=strict
ProtectHome=read-only
ReadWritePaths=/tmp/coolsolve_sessions
PrivateTmp=false
ProtectKernelTunables=true
ProtectControlGroups=true
RestrictSUIDSGID=true

# Environment
Environment=HOME=/home/coolsolve

# Logging
StandardOutput=journal
StandardError=journal
SyslogIdentifier=coolsolve

[Install]
WantedBy=multi-user.target
EOF
```

Enable and start the service:

```bash
sudo systemctl daemon-reload
sudo systemctl enable coolsolve
sudo systemctl start coolsolve

# Verify it's running
sudo systemctl status coolsolve
journalctl -u coolsolve -f    # follow logs
```

Useful commands:

```bash
sudo systemctl restart coolsolve   # restart after recompilation
sudo systemctl stop coolsolve      # stop the server
journalctl -u coolsolve --since today  # today's logs
```

---

## 5. Apache Reverse Proxy Configuration

### 5a. HTTP only (no SSL)

```bash
sudo tee /etc/apache2/sites-available/coolsolve.conf > /dev/null << 'EOF'
<VirtualHost *:80>
    ServerName coolsolve.example.com

    # --- Authentication ---
    <Location />
        AuthType Basic
        AuthName "CoolSolve"
        AuthUserFile /etc/apache2/.htpasswd-coolsolve
        Require valid-user
    </Location>

    # --- Reverse Proxy ---
    ProxyPreserveHost On
    ProxyPass / http://127.0.0.1:8550/
    ProxyPassReverse / http://127.0.0.1:8550/

    # --- SSE support (Server-Sent Events for solve progress) ---
    # Disable buffering for the streaming endpoint
    <Location /api/v1/solve/stream>
        ProxyPass http://127.0.0.1:8550/api/v1/solve/stream
        ProxyPassReverse http://127.0.0.1:8550/api/v1/solve/stream
        # Disable proxy buffering for SSE
        SetEnv proxy-sendchunked 1
        SetEnv proxy-interim-response RFC
    </Location>

    # --- Security Headers ---
    Header always set X-Content-Type-Options "nosniff"
    Header always set X-Frame-Options "DENY"
    Header always set Referrer-Policy "strict-origin-when-cross-origin"

    # --- Logging ---
    ErrorLog ${APACHE_LOG_DIR}/coolsolve_error.log
    CustomLog ${APACHE_LOG_DIR}/coolsolve_access.log combined
</VirtualHost>
EOF
```

### 5b. HTTPS with Let's Encrypt (recommended for production)

Install Certbot:

```bash
sudo apt install -y certbot python3-certbot-apache
```

Create the HTTPS site config:

```bash
sudo tee /etc/apache2/sites-available/coolsolve-ssl.conf > /dev/null << 'EOF'
<VirtualHost *:80>
    ServerName coolsolve.example.com
    # Redirect all HTTP to HTTPS
    RewriteEngine On
    RewriteRule ^(.*)$ https://%{HTTP_HOST}$1 [R=301,L]
</VirtualHost>

<VirtualHost *:443>
    ServerName coolsolve.example.com

    # --- SSL ---
    SSLEngine on
    SSLCertificateFile /etc/letsencrypt/live/coolsolve.example.com/fullchain.pem
    SSLCertificateKeyFile /etc/letsencrypt/live/coolsolve.example.com/privkey.pem

    # --- Authentication ---
    <Location />
        AuthType Basic
        AuthName "CoolSolve"
        AuthUserFile /etc/apache2/.htpasswd-coolsolve
        Require valid-user
    </Location>

    # --- Reverse Proxy ---
    ProxyPreserveHost On
    ProxyPass / http://127.0.0.1:8550/
    ProxyPassReverse / http://127.0.0.1:8550/

    # --- SSE support ---
    <Location /api/v1/solve/stream>
        ProxyPass http://127.0.0.1:8550/api/v1/solve/stream
        ProxyPassReverse http://127.0.0.1:8550/api/v1/solve/stream
        SetEnv proxy-sendchunked 1
        SetEnv proxy-interim-response RFC
    </Location>

    # --- Security Headers ---
    Header always set X-Content-Type-Options "nosniff"
    Header always set X-Frame-Options "DENY"
    Header always set Referrer-Policy "strict-origin-when-cross-origin"
    Header always set Strict-Transport-Security "max-age=63072000; includeSubDomains"

    # --- Logging ---
    ErrorLog ${APACHE_LOG_DIR}/coolsolve_error.log
    CustomLog ${APACHE_LOG_DIR}/coolsolve_access.log combined
</VirtualHost>
EOF
```

Obtain the certificate and enable:

```bash
# Temporarily enable plain HTTP site for Certbot challenge
sudo a2ensite coolsolve.conf
sudo systemctl reload apache2

# Obtain certificate
sudo certbot certonly --apache -d coolsolve.example.com

# Disable plain HTTP, enable SSL site
sudo a2dissite coolsolve.conf
sudo a2ensite coolsolve-ssl.conf
sudo systemctl reload apache2
```

> **Note:** Replace `coolsolve.example.com` with your actual domain name throughout.

---

## 6. User Authentication

### Create the password file and add users

```bash
# Create the file and add the first user (will prompt for password)
sudo htpasswd -c /etc/apache2/.htpasswd-coolsolve admin

# Add additional users (without -c, to avoid overwriting)
sudo htpasswd /etc/apache2/.htpasswd-coolsolve alice
sudo htpasswd /etc/apache2/.htpasswd-coolsolve bob

# Secure the password file
sudo chown root:www-data /etc/apache2/.htpasswd-coolsolve
sudo chmod 640 /etc/apache2/.htpasswd-coolsolve
```

### Manage users

```bash
# Change a user's password
sudo htpasswd /etc/apache2/.htpasswd-coolsolve alice

# Remove a user
sudo htpasswd -D /etc/apache2/.htpasswd-coolsolve bob

# List users
cut -d: -f1 /etc/apache2/.htpasswd-coolsolve
```

> **Important:** Basic Auth sends credentials in base64 (not encrypted). Always use HTTPS in production to protect passwords in transit.

---

## 7. Enable the Site and Test

```bash
# Enable the site
sudo a2ensite coolsolve.conf        # or coolsolve-ssl.conf for HTTPS

# Disable the default site if not needed
# sudo a2dissite 000-default.conf

# Test Apache configuration
sudo apache2ctl configtest

# Reload Apache
sudo systemctl reload apache2
```

Verify everything works:

```bash
# Check CoolSolve is listening
curl -s http://127.0.0.1:8550/api/v1/health
# Expected: {"coolpropReady":true,"status":"ok"}

# Check Apache proxy (with auth)
curl -s -u admin:yourpassword http://coolsolve.example.com/api/v1/health
```

---

## 8. Firewall Configuration

```bash
# Allow HTTP and HTTPS through the firewall
sudo ufw allow 'Apache Full'

# Ensure CoolSolve port is NOT exposed externally (it should only be on 127.0.0.1)
# Verify:
sudo ss -tlnp | grep 8550
# Should show: 127.0.0.1:8550 or 0.0.0.0:8550
```

If CoolSolve binds to `0.0.0.0:8550` (it does by default), block external access to that port:

```bash
sudo ufw deny 8550
```

> **Note:** A future improvement would be to make CoolSolve bind to `127.0.0.1` only (via a `--bind` flag). For now, the firewall rule above prevents external access.

---

## 9. Updating CoolSolve

After recompiling CoolSolve (e.g., after a `git pull` or `svn update`):

```bash
cd /home/coolsolve/CoolSolve/build
sudo -u coolsolve cmake ..
sudo -u coolsolve make -j$(nproc) coolsolve

# Restart the service to load the new binary
sudo systemctl restart coolsolve

# Verify
sudo systemctl status coolsolve
```

If the frontend was also rebuilt:

```bash
cd /home/coolsolve/CoolSolve/gui
sudo -u coolsolve npm run build

cd /home/coolsolve/CoolSolve/build
sudo -u coolsolve cmake ..        # re-embed assets
sudo -u coolsolve make -j$(nproc) coolsolve

sudo systemctl restart coolsolve
```

---

## 10. Running on a Subpath (optional)

If CoolSolve should be available at `https://yourserver.com/coolsolve/` instead of root, the Apache config needs path rewriting. This requires CoolSolve API calls to be aware of the prefix; the current embedded SPA uses relative paths so this should work:

```apache
<Location /coolsolve/>
    AuthType Basic
    AuthName "CoolSolve"
    AuthUserFile /etc/apache2/.htpasswd-coolsolve
    Require valid-user

    ProxyPass http://127.0.0.1:8550/
    ProxyPassReverse http://127.0.0.1:8550/
</Location>
```

---

## 11. Monitoring and Log Rotation

### View logs

```bash
# CoolSolve application logs
journalctl -u coolsolve -f

# Apache access/error logs
tail -f /var/log/apache2/coolsolve_access.log
tail -f /var/log/apache2/coolsolve_error.log
```

### Automatic log rotation

Apache logs are rotated automatically by `logrotate` (installed by default on Ubuntu). CoolSolve logs go through `journald`, which handles its own retention. To configure journald retention:

```bash
sudo tee -a /etc/systemd/journald.conf > /dev/null << 'EOF'
SystemMaxUse=500M
MaxRetentionSec=30day
EOF

sudo systemctl restart systemd-journald
```

---

## 12. Troubleshooting

| Symptom | Check |
|---------|-------|
| 503 Service Unavailable | Apache is up but CoolSolve is not — see below |
| 502 Bad Gateway | `sudo systemctl status coolsolve` — is it running? |
| 403 Forbidden | Check `.htpasswd` file permissions and Apache auth config |
| SSE/progress not streaming | Verify `proxy-sendchunked` env is set in Apache config |
| Slow CoolProp startup | First request after restart triggers warmup (~1-2 s); this is normal |
| Session lost between requests | Ensure Apache forwards cookies (`ProxyPreserveHost On`) |
| Permission denied in logs | Check that `coolsolve` user owns the working directory |
| Service exits with `226/NAMESPACE` | `/tmp/coolsolve_sessions` missing — recreate it (see §3) and ensure `/etc/tmpfiles.d/coolsolve.conf` is installed |

Check all components:

```bash
# 1. Is CoolSolve running?
sudo systemctl status coolsolve

# 2. Can it be reached locally?
curl -s http://127.0.0.1:8550/api/v1/health

# 3. Is Apache running?
sudo systemctl status apache2

# 4. Apache config valid?
sudo apache2ctl configtest

# 5. Firewall rules
sudo ufw status verbose
```

---

## Summary

| Component | Location |
|-----------|----------|
| CoolSolve binary | `/home/coolsolve/CoolSolve/build/coolsolve` |
| Systemd service | `/etc/systemd/system/coolsolve.service` |
| Apache site config | `/etc/apache2/sites-available/coolsolve.conf` |
| Password file | `/etc/apache2/.htpasswd-coolsolve` |
| Application logs | `journalctl -u coolsolve` |
| Apache logs | `/var/log/apache2/coolsolve_{access,error}.log` |
| Session temp data | `/tmp/coolsolve_sessions/` |
| Session dir persistence | `/etc/tmpfiles.d/coolsolve.conf` |
| Examples | `/home/coolsolve/CoolSolve/examples/` |
