#!/bin/bash
# ATGCFLOW ngrok Management Script

BACKEND_DIR="/media/drprabudh/m3/Nextflow-Script/WholeExome/backend"
NGROK_PID_FILE="$BACKEND_DIR/ngrok.pid"
NGROK_LOG_FILE="$BACKEND_DIR/logs/ngrok.log"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to start ngrok
start_ngrok() {
    echo -e "${YELLOW}🚀 Starting ngrok tunnel...${NC}"

    # Check if ngrok is already running
    if [ -f "$NGROK_PID_FILE" ]; then
        NGROK_PID=$(cat "$NGROK_PID_FILE")
        if ps -p $NGROK_PID > /dev/null 2>&1; then
            echo -e "${YELLOW}⚠️  ngrok is already running (PID: $NGROK_PID)${NC}"
            get_ngrok_url
            return
        fi
    fi

    # Kill any existing ngrok processes
    pkill -f ngrok || true
    sleep 1

    # Start ngrok in background
    nohup ngrok http 8000 > "$NGROK_LOG_FILE" 2>&1 &
    NGROK_PID=$!
    echo $NGROK_PID > "$NGROK_PID_FILE"

    echo -e "${GREEN}✅ ngrok started (PID: $NGROK_PID)${NC}"
    echo "   Waiting for ngrok to initialize..."
    sleep 3

    get_ngrok_url
}

# Function to stop ngrok
stop_ngrok() {
    echo -e "${YELLOW}🛑 Stopping ngrok tunnel...${NC}"

    if [ -f "$NGROK_PID_FILE" ]; then
        NGROK_PID=$(cat "$NGROK_PID_FILE")
        if ps -p $NGROK_PID > /dev/null 2>&1; then
            kill $NGROK_PID
            rm "$NGROK_PID_FILE"
            echo -e "${GREEN}✅ ngrok stopped${NC}"
        else
            echo -e "${YELLOW}⚠️  ngrok is not running${NC}"
            rm "$NGROK_PID_FILE"
        fi
    else
        # Try to kill any ngrok process
        pkill -f ngrok && echo -e "${GREEN}✅ ngrok stopped${NC}" || echo -e "${YELLOW}⚠️  No ngrok process found${NC}"
    fi
}

# Function to restart ngrok
restart_ngrok() {
    echo -e "${YELLOW}🔄 Restarting ngrok tunnel...${NC}"
    stop_ngrok
    sleep 1
    start_ngrok
}

# Function to get ngrok URL
get_ngrok_url() {
    NGROK_URL=$(curl -s http://localhost:4040/api/tunnels 2>/dev/null | grep -o '"public_url":"https://[^"]*' | grep -o 'https://[^"]*' | head -1)

    if [ -n "$NGROK_URL" ]; then
        echo -e "${GREEN}✅ ngrok tunnel active!${NC}"
        echo -e "   Public URL: ${GREEN}$NGROK_URL${NC}"
        echo ""
        echo -e "${YELLOW}⚠️  IMPORTANT:${NC} Update your backend/.env file:"
        echo "   Add this URL to CORS_ORIGINS:"
        echo "   $NGROK_URL"
        echo ""
        echo "   Then restart backend:"
        echo "   sudo systemctl restart atgcflow-backend.service"
        echo ""
        echo -e "🔍 ngrok dashboard: ${GREEN}http://localhost:4040${NC}"
    else
        echo -e "${RED}❌ Could not get ngrok URL${NC}"
        echo "   Check manually: http://localhost:4040"
        echo "   Or view logs: tail -f $NGROK_LOG_FILE"
    fi
}

# Function to show status
show_status() {
    echo -e "${YELLOW}📊 ngrok Status${NC}"
    echo "================="
    echo ""

    if [ -f "$NGROK_PID_FILE" ]; then
        NGROK_PID=$(cat "$NGROK_PID_FILE")
        if ps -p $NGROK_PID > /dev/null 2>&1; then
            echo -e "Status: ${GREEN}Running${NC} (PID: $NGROK_PID)"
        else
            echo -e "Status: ${RED}Not Running${NC} (stale PID file)"
            rm "$NGROK_PID_FILE"
        fi
    else
        # Check if ngrok is running without PID file
        NGROK_RUNNING=$(pgrep -f "ngrok http")
        if [ -n "$NGROK_RUNNING" ]; then
            echo -e "Status: ${GREEN}Running${NC} (PID: $NGROK_RUNNING - no PID file)"
        else
            echo -e "Status: ${RED}Not Running${NC}"
        fi
    fi

    echo ""
    get_ngrok_url
}

# Main script
case "$1" in
    start)
        start_ngrok
        ;;
    stop)
        stop_ngrok
        ;;
    restart)
        restart_ngrok
        ;;
    status)
        show_status
        ;;
    url)
        get_ngrok_url
        ;;
    *)
        echo "Usage: $0 {start|stop|restart|status|url}"
        echo ""
        echo "Commands:"
        echo "  start   - Start ngrok tunnel"
        echo "  stop    - Stop ngrok tunnel"
        echo "  restart - Restart ngrok tunnel"
        echo "  status  - Show ngrok status"
        echo "  url     - Get current ngrok URL"
        echo ""
        exit 1
        ;;
esac

exit 0
