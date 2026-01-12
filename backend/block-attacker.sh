#!/bin/bash
# Emergency script to block malicious IP

ATTACKER_IP="192.168.1.27"

echo "🚨 BLOCKING MALICIOUS IP: $ATTACKER_IP"
echo ""

# Block using iptables
echo "Adding iptables rule..."
sudo iptables -I INPUT -s $ATTACKER_IP -j DROP
sudo iptables -I OUTPUT -d $ATTACKER_IP -j DROP

echo "✅ IP blocked in iptables"
echo ""

# Make persistent (Ubuntu/Debian)
if command -v iptables-save &> /dev/null; then
    echo "Making iptables rules persistent..."
    sudo netfilter-persistent save 2>/dev/null || sudo iptables-save | sudo tee /etc/iptables/rules.v4 > /dev/null
    echo "✅ Rules saved"
fi

echo ""
echo "📋 Current iptables rules:"
sudo iptables -L INPUT -v -n | grep $ATTACKER_IP

echo ""
echo "✅ Attacker IP blocked successfully!"
echo ""
echo "To unblock later (if needed):"
echo "  sudo iptables -D INPUT -s $ATTACKER_IP -j DROP"
echo "  sudo iptables -D OUTPUT -d $ATTACKER_IP -j DROP"
