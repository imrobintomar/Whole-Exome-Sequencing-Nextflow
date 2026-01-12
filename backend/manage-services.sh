#!/bin/bash
# Unified service management script for ATGCFLOW Backend + ngrok

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_header() {
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}  $1${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

get_ngrok_url() {
    curl -s http://localhost:4040/api/tunnels 2>/dev/null | grep -o '"public_url":"https://[^"]*' | grep -o 'https://[^"]*' | head -1
}

check_service_status() {
    local service=$1
    if systemctl is-active --quiet $service; then
        print_success "$service is running"
        return 0
    else
        print_error "$service is not running"
        return 1
    fi
}

status_command() {
    print_header "ATGCFLOW Services Status"

    echo "Backend Service:"
    check_service_status "atgcflow-backend.service"
    if systemctl is-enabled --quiet atgcflow-backend.service 2>/dev/null; then
        print_info "Auto-start: enabled"
    else
        print_warning "Auto-start: disabled"
    fi
    echo ""

    echo "ngrok Service:"
    check_service_status "atgcflow-ngrok.service"
    if systemctl is-enabled --quiet atgcflow-ngrok.service 2>/dev/null; then
        print_info "Auto-start: enabled"
    else
        print_warning "Auto-start: disabled"
    fi
    echo ""

    # Check if backend is responding
    if curl -s http://localhost:8000/ > /dev/null 2>&1; then
        print_success "Backend API responding at http://localhost:8000"
    else
        print_warning "Backend API not responding"
    fi

    # Get ngrok URL
    NGROK_URL=$(get_ngrok_url)
    if [ -n "$NGROK_URL" ]; then
        print_success "ngrok tunnel active: $NGROK_URL"
        echo "$NGROK_URL" > "$SCRIPT_DIR/ngrok-url.txt"
    else
        print_warning "ngrok tunnel not found (check http://localhost:4040)"
    fi

    echo ""
}

start_command() {
    print_header "Starting ATGCFLOW Services"

    echo "Starting backend service..."
    sudo systemctl start atgcflow-backend.service
    sleep 2
    check_service_status "atgcflow-backend.service"
    echo ""

    echo "Starting ngrok service..."
    sudo systemctl start atgcflow-ngrok.service
    sleep 3
    check_service_status "atgcflow-ngrok.service"
    echo ""

    NGROK_URL=$(get_ngrok_url)
    if [ -n "$NGROK_URL" ]; then
        print_success "Public URL: $NGROK_URL"
        echo "$NGROK_URL" > "$SCRIPT_DIR/ngrok-url.txt"
    fi
}

stop_command() {
    print_header "Stopping ATGCFLOW Services"

    echo "Stopping ngrok service..."
    sudo systemctl stop atgcflow-ngrok.service
    print_success "ngrok stopped"
    echo ""

    echo "Stopping backend service..."
    sudo systemctl stop atgcflow-backend.service
    print_success "Backend stopped"
    echo ""
}

restart_command() {
    print_header "Restarting ATGCFLOW Services"

    echo "Restarting backend service..."
    sudo systemctl restart atgcflow-backend.service
    sleep 2
    check_service_status "atgcflow-backend.service"
    echo ""

    echo "Restarting ngrok service..."
    sudo systemctl restart atgcflow-ngrok.service
    sleep 3
    check_service_status "atgcflow-ngrok.service"
    echo ""

    NGROK_URL=$(get_ngrok_url)
    if [ -n "$NGROK_URL" ]; then
        print_success "New public URL: $NGROK_URL"
        echo "$NGROK_URL" > "$SCRIPT_DIR/ngrok-url.txt"
        print_warning "ngrok URL has changed! Update your frontend configuration."
    fi
}

logs_command() {
    local service=$1

    case $service in
        backend)
            print_info "Following backend logs (Ctrl+C to stop)..."
            sudo journalctl -u atgcflow-backend.service -f
            ;;
        ngrok)
            print_info "Following ngrok logs (Ctrl+C to stop)..."
            sudo journalctl -u atgcflow-ngrok.service -f
            ;;
        all)
            print_info "Following all logs (Ctrl+C to stop)..."
            sudo journalctl -u atgcflow-backend.service -u atgcflow-ngrok.service -f
            ;;
        *)
            print_error "Invalid service. Use: backend, ngrok, or all"
            exit 1
            ;;
    esac
}

url_command() {
    NGROK_URL=$(get_ngrok_url)
    if [ -n "$NGROK_URL" ]; then
        echo "$NGROK_URL"
        # Also copy to clipboard if xclip is available
        if command -v xclip &> /dev/null; then
            echo "$NGROK_URL" | xclip -selection clipboard
            print_success "URL copied to clipboard"
        fi
    else
        print_error "ngrok tunnel not active"
        print_info "Check: http://localhost:4040"
        exit 1
    fi
}

enable_command() {
    print_header "Enabling Auto-start on Boot"

    sudo systemctl enable atgcflow-backend.service
    print_success "Backend auto-start enabled"

    sudo systemctl enable atgcflow-ngrok.service
    print_success "ngrok auto-start enabled"
    echo ""
}

disable_command() {
    print_header "Disabling Auto-start on Boot"

    sudo systemctl disable atgcflow-backend.service
    print_warning "Backend auto-start disabled"

    sudo systemctl disable atgcflow-ngrok.service
    print_warning "ngrok auto-start disabled"
    echo ""
}

usage() {
    cat << EOF
ATGCFLOW Service Manager

Usage: $0 <command> [options]

Commands:
    status              Show status of all services
    start               Start all services
    stop                Stop all services
    restart             Restart all services
    logs <service>      Follow logs (service: backend, ngrok, or all)
    url                 Get current ngrok public URL
    enable              Enable auto-start on boot
    disable             Disable auto-start on boot

Examples:
    $0 status                # Check service status
    $0 start                 # Start all services
    $0 logs backend          # View backend logs
    $0 logs all              # View all logs
    $0 url                   # Get ngrok URL
    $0 restart               # Restart everything

EOF
}

# Main command router
case "${1:-}" in
    status)
        status_command
        ;;
    start)
        start_command
        ;;
    stop)
        stop_command
        ;;
    restart)
        restart_command
        ;;
    logs)
        logs_command "${2:-all}"
        ;;
    url)
        url_command
        ;;
    enable)
        enable_command
        ;;
    disable)
        disable_command
        ;;
    help|--help|-h)
        usage
        ;;
    *)
        print_error "Invalid command: ${1:-}"
        echo ""
        usage
        exit 1
        ;;
esac
